/*
 * This code is provided solely for the personal and private use of students
 * taking the CSC369H course at the University of Toronto. Copying for purposes
 * other than this use is expressly prohibited. All forms of distribution of
 * this code, including but not limited to public repositories on GitHub,
 * GitLab, Bitbucket, or any other online platform, whether as given or with
 * any changes, are expressly prohibited.
 *
 * Authors: Alexey Khrabrov, Karen Reid
 *
 * All of the files in this directory and all subdirectories are:
 * Copyright (c) 2019 Karen Reid
 */

/**
 * CSC369 Assignment 1 - a1fs driver implementation.
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>

// Using 2.9.x FUSE API
#define FUSE_USE_VERSION 29
#include <fuse.h>

#include "a1fs.h"
#include "fs_ctx.h"
#include "options.h"
#include "map.h"

#include <time.h>
#include <string.h>

static fs_ctx *get_fs(void);
//NOTE: All path arguments are absolute paths within the a1fs file system and
// start with a '/' that corresponds to the a1fs root directory.
//
// For example, if a1fs is mounted at "~/my_csc369_repo/a1b/mnt/", the path to a
// file at "~/my_csc369_repo/a1b/mnt/dir/file" (as seen by the OS) will be
// passed to FUSE callbacks as "/dir/file".
//
// Paths to directories (except for the root directory - "/") do not end in a
// trailing '/'. For example, "~/my_csc369_repo/a1b/mnt/dir/" will be passed to
// FUSE callbacks as "/dir".

/**
 * Set the bitmap bit to 0 according to num.
 * 
 * @param ibmp  1 if the bitmap to be set is inode bitmap, 0 if the bitmap to be set is block bitmap.
 * @param num   the position of the bit to be set, starting from 0.
 */
void set_bitmap_to_0(int ibmp, int num){
	fs_ctx *fs = get_fs();
	struct a1fs_superblock *sb = fs->image;
	unsigned char *start;
	if(ibmp) start = fs->image + sb->inode_bitmap * A1FS_BLOCK_SIZE;
	else start = fs->image + sb->block_bitmap * A1FS_BLOCK_SIZE;	
	int byte = num / 8;
	int bit = num % 8;
	start[byte] = start[byte] & ~(1 << bit);
	if(ibmp) sb->free_inode_count += 1;
	else sb->free_block_count += 1;
}


/**
 * Find the first 0 in inode/block bitmap and set to 1.
 * 
 * @param ibmp  1 if the bitmap to be set is inode bitmap, 0 if the bitmap to be set is block bitmap.
 * @return      the number indicating the position of the bit just set, 0 if failed to set any bit.
 */
int set_first_available(int ibmp){
	fs_ctx *fs = get_fs();
	struct a1fs_superblock *sb = fs->image;
	int num_bytes;
	unsigned char *start;
	int remainder;
	int flag = 0; 
	if(ibmp){
		start = fs->image + sb->inode_bitmap * A1FS_BLOCK_SIZE;
		num_bytes = sb->inode_count / 8;
		remainder = sb->inode_count % 8;
	} 
	else{
		start = fs->image + sb->block_bitmap * A1FS_BLOCK_SIZE;
		num_bytes = sb->block_count / 8;
		remainder = sb->block_count % 8;	
	}
	for(int byte = 0; byte < num_bytes; byte++){
		if(start[byte] != 255){
			for(int bit = 0; bit < 8; bit++){
				if(!(start[byte] & (1 << bit))){
					start[byte] = start[byte] | (1 << bit);
					if(ibmp) sb->free_inode_count -= 1;
					else sb->free_block_count -= 1;
					return byte * 8 + bit;
				} 
			}
		}
		flag = 1; //indicates the for loop has been run.
	}
	if(flag) num_bytes += 1; //accomodate the special case: no enough bits to form an entire byte.
	if(remainder > 0){ //excess bits
		for(int i = 0; i < remainder; i++){
			if(!(start[num_bytes] & (1 << i))){
				start[num_bytes] = start[num_bytes] | (1 << i);
				if(ibmp) sb->free_inode_count -= 1;
				else sb->free_block_count -= 1;
				return (num_bytes) * 8 + i;
			}
		}
	}
	return 0;
}


/**
 * Find the extent of the file containing the byte indicated by offset.
 * 
 * @param inode      pointer to the inode of the file.
 * @param offset     the nth byte to be found, starting from 0.
 * @param res_ext    pointer to the integer indicating the nth extent from 0, which is the ext storing the target.
 * @param res_blk    pointer to the integer indicating the nth block form 0, which is the blk storing the target.
 * @param remainder  pointer to the position of the (offset)th byte in the target extent.
 * @return           0 on success, 1 on failure.
 */
int find_ext_by_offset(struct a1fs_inode *inode, int offset, unsigned int *res_ext, int *res_blk, int *remainder){
	fs_ctx *fs = get_fs();
	unsigned int remain = offset;
	struct a1fs_extent *ext = inode->extents;
	unsigned int count = 1; //number of ext, starting from 1
	while(remain > A1FS_BLOCK_SIZE * ext->count){ //find the ext containing (offset)th byte
		remain -= A1FS_BLOCK_SIZE * ext->count;
		if(count < 10 || (count > 10 && count < 512)) ext = ext + 1; //can safely move to the next ext
		if(count == 10) ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE; //jump to extent_blk
		count += 1;
		if(count > inode->extent_count) return 1;
	}
	*remainder = remain % A1FS_BLOCK_SIZE;
	*res_blk = remain / A1FS_BLOCK_SIZE;
	*res_ext = count - 1; //starts from 0
	return 0;
}


/**
 * Free all blocks in the extent and update bitmap block.
 * 
 * @param ext  pointer to the extent.
 */
void reset_blocks_in_ext(struct a1fs_extent *ext){
	fs_ctx *fs = get_fs();
	for(unsigned int i = ext->start; i < ext->start + ext->count; i++){
		set_bitmap_to_0(0, i);
	}
	memset(fs->image + ext->start * A1FS_BLOCK_SIZE, 0, ext->count * A1FS_BLOCK_SIZE);
	ext->start = 0;
	ext->count = 0;
}


/**
 * Free all blocks in the inode starting from start_ext. Free the block storing the extra extents as well.
 * 
 * @param  inode      pointer to the inode.
 * @param  start_ext  the index of the first extent to be reset.
 */
void reset_blocks_in_ino(struct a1fs_inode *inode, int start_ext){
	fs_ctx *fs = get_fs();
	if(start_ext < 10){ //need to reset some extents in the array
		for(unsigned int i = start_ext; i < 10; i++){
			if(i < inode->extent_count){
				reset_blocks_in_ext(inode->extents + i);
			}
		}
		if(inode->extent_count > 10){ //then reset extents from the extra block
			struct a1fs_extent *ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE;
			for(unsigned int j = 0; j < inode->extent_count - 10; j++){
				reset_blocks_in_ext(ext + j);
			}
		}
	}
	else{ //start from the extra block
		if(inode->extent_count > 10){
			struct a1fs_extent *ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE + (start_ext-10) * sizeof(struct a1fs_extent);
			for(unsigned int j = start_ext; j < inode->extent_count - 10; j++){
				reset_blocks_in_ext(ext + j);
			}
		}
	}
	if(inode->extent_blk != 0 && start_ext <= 10){
		//reset the block storing extra extents
		struct a1fs_superblock *sb = fs->image;
		unsigned char *block_bitmap = fs->image + sb->block_bitmap;
		int byte = inode->extent_blk / 8;
		int bit = inode->extent_blk % 8;
		block_bitmap[byte] = block_bitmap[byte] & ~(1 << bit);
		sb->free_block_count += 1;
		memset(fs->image + inode->extent_blk * A1FS_BLOCK_SIZE, 0, A1FS_BLOCK_SIZE);
	}
	inode->extent_count = start_ext;
}


/**
 * Shrink the file to the given size.
 * 
 * @param inode  pointer to the inode of the file.
 * @param size   the size to be shrinked to.
 */
void shrink_file(struct a1fs_inode *inode, int size){
	fs_ctx *fs = get_fs();
	a1fs_superblock *sb = fs->image;
	struct a1fs_extent *first_ext = inode->extents;
	unsigned int count = 0;
	int blocks = 0;
	int remainder = 0;
	find_ext_by_offset(inode, size, &count, &blocks, &remainder);
	if(count <= 9) first_ext = inode->extents + count;
	else first_ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE + (count-10) * sizeof(struct a1fs_extent);
	count += 1;
	int bit;
	//shrink the ext just found
	memset(fs->image + first_ext->start * A1FS_BLOCK_SIZE + blocks * A1FS_BLOCK_SIZE + remainder, 0, first_ext->count * A1FS_BLOCK_SIZE - blocks * A1FS_BLOCK_SIZE - remainder);
	unsigned int start = first_ext->start + blocks;
	if(remainder > 0) start += 1;
	for(unsigned int i = start; i < first_ext->start + first_ext->count; i++){
		set_bitmap_to_0(0, i);
	}
	if(remainder > 0) blocks += 1;
	first_ext->count = blocks;
	reset_blocks_in_ino(inode, count); //clear all ext after first_ext
	if(blocks == 0 && remainder == 0){ //first_ext is also cleared
		first_ext->start = 0;
		first_ext->count = 0;
		inode->extent_count -= 1;
		if(count == 11){ //first_ext has index 10: reset the extent_blk(not done in reset_blocks_in_ino)
			unsigned char *block_bitmap = fs->image + sb->block_bitmap;
			int byte = inode->extent_blk / 8;
			bit = inode->extent_blk % 8;
			block_bitmap[byte] = block_bitmap[byte] & ~(1 << bit);
			sb->free_block_count += 1;
			memset(fs->image + inode->extent_blk * A1FS_BLOCK_SIZE, 0, A1FS_BLOCK_SIZE);
		}
	}
}


/**
 * Expand the file to the given size.
 * 
 * @param inode  pointer to the inode of the file.
 * @param size   the size to be expanded to.
 * @return       0 on success, 1 on failure.
 */
int expand_file(struct a1fs_inode *inode, int size){
	fs_ctx *fs = get_fs();
	a1fs_superblock *sb = fs->image;
	int remain = size;
	struct a1fs_extent *last_ext;
	//find the last extent of the file
	if(inode->extent_count <= 10) last_ext = (inode->extents) + inode->extent_count - 1;
	else last_ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE + (inode->extent_count-1-10) * sizeof(struct a1fs_extent);
	int byte = inode->size % A1FS_BLOCK_SIZE;
	int bit;
	//expand the last extent
	if(byte != 0 || inode->size == 0){
		if(remain <= A1FS_BLOCK_SIZE - byte){ //size to expand <= space left
			memset(fs->image + (last_ext->start + last_ext->count - 1) * A1FS_BLOCK_SIZE + byte, 0, remain);
			inode->size += remain;
			return 0;
		}
		else{ //size to expand > space left
			memset(fs->image + (last_ext->start + last_ext->count - 1) * A1FS_BLOCK_SIZE + byte, 0, A1FS_BLOCK_SIZE - byte);
			remain -= A1FS_BLOCK_SIZE - byte;
			inode->size += A1FS_BLOCK_SIZE - byte;
		}
	}
	if(sb->free_block_count == 0) return 1;
	unsigned char *bmp = fs->image + sb->block_bitmap * A1FS_BLOCK_SIZE;
	int last_expanded = 0; //0 indicates the last ext is not expanded, 1 otherwise
	int start;
	struct a1fs_extent *next;
	while(remain > 0){
		if(!last_expanded){ //check if the last ext is expanded, if not, expand
			byte = (last_ext->start + last_ext->count) / 8;
			bit = (last_ext->start + last_ext->count) % 8;
			if(bmp[byte] & (1 << bit)){ //the next block is taken: cannot expand
				last_expanded = 1; //mark expanded
				continue;
			}
			//can expand: take the next block, do not mark expanded, continue checking next time
			last_ext->count += 1;
			bmp[byte] = bmp[byte] | (1 << bit);
			sb->free_block_count -= 1;
			if(remain <= A1FS_BLOCK_SIZE){
				memset(fs->image + (last_ext->start + last_ext->count - 1) * A1FS_BLOCK_SIZE, 0, remain);
				inode->size += remain;
			}
			else{
				memset(fs->image + (last_ext->start + last_ext->count - 1) * A1FS_BLOCK_SIZE, 0, A1FS_BLOCK_SIZE);
				inode->size += A1FS_BLOCK_SIZE;
			}
			remain -= A1FS_BLOCK_SIZE;
		}
		else{ //add new extent
			//less than 10 or less than 512: can directly add one extent
			start = set_first_available(0);
			if(start == 0) return 1;
			if(inode->extent_count < 10){
				next = inode->extents + inode->extent_count;
			}
			else if(inode->extent_count > 10 && inode->extent_count < 512){
				next = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE + (inode->extent_count-10-1) * sizeof(struct a1fs_extent);
			}
			//exactly 10: allocate a new data block for the new extent
			else if(inode->extent_count == 10){
				if((inode->extent_blk = set_first_available(0)) == 0){
					//undo
					byte = start / 8;
					bit = start % 8;
					unsigned char *bmp_start = fs->image + sb->block_bitmap * A1FS_BLOCK_SIZE;
					bmp_start[byte] = bmp_start[byte] & ~(1 << bit);
					sb->free_block_count += 1;
					return 1;
				}
				next = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE;
			}
			else return 1; //exactly 512: cannot fit
			next->start = start;
			next->count = 1;
			if(remain <= A1FS_BLOCK_SIZE){
				memset(fs->image + start * A1FS_BLOCK_SIZE, 0, remain);
				inode->size += remain;
			}
			else{
				memset(fs->image + start * A1FS_BLOCK_SIZE, 0, A1FS_BLOCK_SIZE);
				inode->size += A1FS_BLOCK_SIZE;
			}
			last_ext = next; //the newly added ext is now the last
			inode->extent_count += 1;
			last_expanded = 0; //the new last ext is not expanded
			remain -= A1FS_BLOCK_SIZE;
		}
	}
	return 0;
}


/**
 * Delete a dentry in extent if exists.
 * 
 * @param name  the name of the dentry to be deleted.
 * @param ext   pointer to the extent.
 * @return      the inode number of the dentry file, 0 if such dentry is not found.
 */
int delete_dentry_in_ext(char *name, struct a1fs_extent *ext){
	fs_ctx *fs = get_fs();
	struct a1fs_dentry *first = fs->image + ext->start * A1FS_BLOCK_SIZE;
	for(unsigned int i = 0; i < A1FS_BLOCK_SIZE * ext->count / sizeof(struct a1fs_dentry); i++){
		if(strcmp((first + i)->name, name) == 0){ //check if name matchs
			int ino_num = (first + i)->ino;
			memset(first + i, 0, sizeof(struct a1fs_dentry));
			return ino_num;
		}
	}
	return 0;
}


/**
 * Delete a dentry in inode if exists.
 * 
 * @param inode  pointer to the inode.
 * @param name   the name of the dentry to be deleted.
 * @return       the inode number of the dentry file, 0 if such dentry is not found.
 */
int delete_dentry_in_ino(struct a1fs_inode *inode, char *name){
	fs_ctx *fs = get_fs();
	int res;
	for(unsigned int i = 0; i < 10; i++){
		if(i < inode->extent_count){
			if((res = delete_dentry_in_ext(name, inode->extents + i)) != 0){
				inode->size -= sizeof(struct a1fs_dentry);
				return res;
			}
		}
	}
	if(inode->extent_count > 10){
		struct a1fs_extent *ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE;
		for(unsigned int j = 0; j < inode->extent_count - 10; j++){
			if((res = delete_dentry_in_ext(name, ext + j)) != 0){
				inode->size -= sizeof(struct a1fs_dentry);
				return res;
			}
		}
	}
	return 0;
}


/**
 * Find free entry space in the extent, if cannot find, try expanding the extent.
 * 
 * @param entry  pointer to pointer to the free entry space found.
 * @param ext    pointer to the extent.
 * @return       0 on success, 1 on failure.
 */
int find_dspace_in_ext(struct a1fs_dentry **entry, struct a1fs_extent *ext){
	fs_ctx *fs = get_fs();
	struct a1fs_superblock *sb = fs->image;
	//find available space in the extent or find dentry with target name
	struct a1fs_dentry *first = fs->image + ext->start * A1FS_BLOCK_SIZE;
	for(unsigned int i = 0; i < A1FS_BLOCK_SIZE * ext->count / sizeof(struct a1fs_dentry); i++){
		if((first + i)->ino == 0 && strlen((first + i)->name) == 0){
			*entry = first + i;
			return 0;
		}
	}
	//expand the extent to include the next block
	unsigned char *block_bitmap = fs->image + sb->block_bitmap;
	int target = ext->start + ext->count;
	int byte = target / 8;
	int bit = target % 8;
	if(!(block_bitmap[byte] & (1 << bit))){
		ext->count += 1;
		block_bitmap[byte] = block_bitmap[byte] | (1 << bit);
		*entry = fs->image + A1FS_BLOCK_SIZE * (ext->start + ext->count);
		return 0;
	}
	return 1;
}


/**
 * Find free entry space in the inode.
 * 
 * @param inode  pointer to the inode.
 * @param entry  pointer to pointer to the free entry space found.
 * @return       0 on success, 1 on failure.
 */
int find_despace_in_ino(struct a1fs_inode *inode, struct a1fs_dentry **entry){
	fs_ctx *fs = get_fs();
	for(unsigned int i = 0; i < 10; i++){
		if(i < inode->extent_count){
			if(find_dspace_in_ext(entry, inode->extents + i) == 0) return 0;
		}
	}
	if(inode->extent_count > 10){
		struct a1fs_extent *ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE;
		for(unsigned int j = 0; j < inode->extent_count - 10; j++){
			if(find_dspace_in_ext(entry, ext + j) == 0) return 0;
		}
	}
	return 1;
}


/**
 * Find the dentry with the specified name in the extent.
 * 
 * @param extent pointer to the extent.
 * @param name   the name of the target dentry.
 * @return       the pointer to the dentry found. 
 */
struct a1fs_dentry* find_entry_in_ext(struct a1fs_extent extent, char *name){
	fs_ctx *fs = get_fs();
	struct a1fs_dentry *target = fs->image + extent.start * A1FS_BLOCK_SIZE;
	for(unsigned int i = 0; i < extent.count; i++){
		target = fs->image + extent.start * A1FS_BLOCK_SIZE + i * A1FS_BLOCK_SIZE;
		for(unsigned int j = 0; j < A1FS_BLOCK_SIZE / sizeof(struct a1fs_dentry); j++){
			if(strcmp((target + j)->name, name) == 0){
				return target + j;
			}
		}
	}
	return NULL;
}


/**
 * Find the dentry with the specified name in the inode.
 * 
 * @param inode  pointer to the inode.
 * @param name   the name of the target dentry.
 * @return       the pointer to the dentry found. 
 */
struct a1fs_dentry* find_entry_in_ino(struct a1fs_inode *inode, char *name){
	fs_ctx *fs = get_fs();
	struct a1fs_dentry* target = NULL;
	for(unsigned int i = 0; i < 10; i++){
		if(i < inode->extent_count){
			target = find_entry_in_ext(inode->extents[i], name);
			if(target != NULL) {
				return target;
			}
		}
	}
	if(inode->extent_count > 10){
		struct a1fs_extent *ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE;
		for(unsigned int j = 0; j < inode->extent_count - 10; j++){
			target = find_entry_in_ext(ext[j], name);
			if(target != NULL) return target;
		}
	}
	return NULL;
}


/**
 * Find the inode of the file indicated by the path given.
 * 
 * @param path   the path of the file.
 * @param inode  the pointer to pointer to the target inode.
 * @return       0 on success, 1 if path is not a directory, 2 if such entry is not found.
 */
int find_inode(const char *path, struct a1fs_inode **inode){
	fs_ctx *fs = get_fs();
	struct a1fs_superblock *sb = fs->image;
	char delim[2] = "/";
	char path_name[PATH_MAX];
	strncpy(path_name, path, PATH_MAX);
	char *name = strtok(path_name, delim); //separate the path by '/'
	struct a1fs_inode *curr = fs->image + sb->inode_table * A1FS_BLOCK_SIZE;
	if(strcmp(path, delim) == 0){ //check if the file to be found is the root
		*inode = curr;
		return 0;
	}
	struct a1fs_dentry *target;
	while(name != NULL){ 
		if((target = find_entry_in_ino(curr, name)) == NULL) return 2; //entry not found
		curr = fs->image + sb->inode_table * A1FS_BLOCK_SIZE + target->ino * sizeof(struct a1fs_inode);
		name = strtok(NULL, delim); //go to next level
		if(!S_ISDIR(curr->mode) && name != NULL) return 1; //path is not a dir
	}
	*inode = curr;
	return 0;
}


/**
 * Add a dentry to the inode.
 * 
 * @param inode  pointer to the inode.
 * @param entry  pointer to pointer to the dentry to be added.
 * @return       0 on success, 1 on failure.
 */
int add_dentry(struct a1fs_inode *inode, struct a1fs_dentry **entry){
	fs_ctx *fs = get_fs();
	struct a1fs_extent *ext;
	//try to insert dentry into existing extents or extend existing extents
	if(find_despace_in_ino(inode, entry) == 1){ //all existing extents cannot be extended to include the new dentry.
		//less than 10 or less than 512: can directly add one extent
		if(inode->extent_count < 10 || (inode->extent_count > 10 && inode->extent_count < 512)){
			inode->extent_count += 1;
			if(inode->extent_count <= 10) ext = inode->extents + inode->extent_count - 1;
			else ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE + (inode->extent_count-10-1) * sizeof(struct a1fs_extent);
			find_dspace_in_ext(entry, ext);
		}
		//exactly 10: allocate a new data block for the new extent
		else if(inode->extent_count == 10){
			if((inode->extent_blk = set_first_available(0)) == 0) return 1;
			inode->extent_count += 1;
			ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE;
			find_dspace_in_ext(entry, ext);
		}
		else return 1; //exactly 512: cannot fit
	}
	return 0;
}


/**
 * Call filler on all entries in the extent.
 * 
 * @param extent  pointer to the extent.
 * @param buf     buffer that receives the result.
 * @param filler  function that needs to be called for each directory entry.
 * @return        0 on success, 1 on failure.
 */
int filler_entry_in_ext(struct a1fs_extent extent, void *buf, fuse_fill_dir_t filler){
	fs_ctx *fs = get_fs();
	struct a1fs_dentry *target = fs->image + extent.start * A1FS_BLOCK_SIZE;
	for(unsigned int i = 0; i < extent.count; i++){
		target = fs->image + extent.start * A1FS_BLOCK_SIZE + i * A1FS_BLOCK_SIZE;
		for(unsigned int j = 0; j < A1FS_BLOCK_SIZE / sizeof(struct a1fs_dentry); j++){
			if((target + j)->name[0] != 0){
				if(filler(buf, target[j].name, NULL, 0)) return 1;
			}
		}
	}
	return 0;
}


/**
 * Initialize the file system.
 *
 * Called when the file system is mounted. NOTE: we are not using the FUSE
 * init() callback since it doesn't support returning errors. This function must
 * be called explicitly before fuse_main().
 *
 * @param fs    file system context to initialize.
 * @param opts  command line options.
 * @return      true on success; false on failure.
 */
static bool a1fs_init(fs_ctx *fs, a1fs_opts *opts)
{
	// Nothing to initialize if only printing help or version
	if (opts->help || opts->version) return true;

	size_t size;
	void *image = map_file(opts->img_path, A1FS_BLOCK_SIZE, &size);
	if (!image) return false;

	return fs_ctx_init(fs, image, size, opts);
}

/**
 * Cleanup the file system.
 *
 * Called when the file system is unmounted. Must cleanup all the resources
 * created in a1fs_init().
 */
static void a1fs_destroy(void *ctx)
{
	fs_ctx *fs = (fs_ctx*)ctx;
	if (fs->image) {
		if (fs->opts->sync && (msync(fs->image, fs->size, MS_SYNC) < 0)) {
			perror("msync");
		}
		munmap(fs->image, fs->size);
		fs_ctx_destroy(fs);
	}
}

/** Get file system context. */
static fs_ctx *get_fs(void)
{
	return (fs_ctx*)fuse_get_context()->private_data;
}


/**
 * Get file system statistics.
 *
 * Implements the statvfs() system call. See "man 2 statvfs" for details.
 * The f_bfree and f_bavail fields should be set to the same value.
 * The f_ffree and f_favail fields should be set to the same value.
 * The following fields can be ignored: f_fsid, f_flag.
 * All remaining fields are required.
 *
 * @param path  path to any file in the file system. Can be ignored.
 * @param st    pointer to the struct statvfs that receives the result.
 * @return      0 on success; -errno on error.
 */
static int a1fs_statfs(const char *path, struct statvfs *st)
{	
	// add error check
	if(strlen(path) >= PATH_MAX) return -ENAMETOOLONG;
	fs_ctx *fs = get_fs();

	memset(st, 0, sizeof(*st));
	st->f_bsize   = A1FS_BLOCK_SIZE;
	st->f_frsize  = A1FS_BLOCK_SIZE;
	//TODO: fill in the rest of required fields based on the information stored
	// in the superblock
	struct a1fs_superblock *sb = (struct a1fs_superblock *)fs->image;
	st->f_blocks = sb->block_count;
	st->f_bfree = sb->free_block_count;
	st->f_bavail = sb->free_block_count;
	st->f_files = sb->inode_count;
	st->f_ffree = sb->free_inode_count;
	st->f_favail = sb->free_inode_count;
	st->f_namemax = A1FS_NAME_MAX;
	return 0;
}


/**
 * Get file or directory attributes.
 *
 * Implements the stat() system call. See "man 2 stat" for details.
 * The following fields can be ignored: st_dev, st_ino, st_uid, st_gid, st_rdev,
 *                                      st_blksize, st_atim, st_ctim.
 * All remaining fields are required.
 *
 * NOTE: the st_blocks field is measured in 512-byte units (disk sectors).
 *
 * Errors:
 *   ENAMETOOLONG  the path or one of its components is too long.
 *   ENOENT        a component of the path does not exist.
 *   ENOTDIR       a component of the path prefix is not a directory.
 *
 * @param path  path to a file or directory.
 * @param st    pointer to the struct stat that receives the result.
 * @return      0 on success; -errno on error;
 */
static int a1fs_getattr(const char *path, struct stat *st)
{
	if (strlen(path) >= A1FS_PATH_MAX) return -ENAMETOOLONG;
	fs_ctx *fs = get_fs();

	memset(st, 0, sizeof(*st));

	//TODO: lookup the inode for given path and, if it exists, fill in the
	// required fields based on the information stored in the inode
	char *last_delim = strrchr(path, '/'); //last occurrence
	char parent[NAME_MAX];
	strncpy(parent, path, last_delim - path + 1);
	parent[last_delim - path + 1] = '\0'; //get parent path
	struct a1fs_inode *pinode;
	int pres = find_inode(parent, &pinode);
	if(pres == 1) return -ENOTDIR;
	if(pres == 2) return -ENOENT;
	struct a1fs_inode *cinode;
	int cres = find_inode(path, &cinode);
	if(cres == 2) return -ENOENT;
	st->st_mode = cinode->mode;
	st->st_nlink = cinode->links;
	st->st_size = cinode->size;
	int res = 0;
	//count blocks in all extents
	for(unsigned int i = 0; i < 10; i++){
		if(i < cinode->extent_count){
			res += cinode->extents[i].count;
		}
	}
	if(cinode->extent_count > 10){
		struct a1fs_extent *ext = fs->image + cinode->extent_blk * A1FS_BLOCK_SIZE;
		for(unsigned int j = 0; j < cinode->extent_count - 10; j++){
			res += ext[j].count;
		}
	}
	st->st_blocks = res;
	st->st_mtim = cinode->mtime;
	return 0;
}


/**
 * Read a directory.
 *
 * Implements the readdir() system call. Should call filler() for each directory
 * entry. See fuse.h in libfuse source code for details.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a directory.
 *
 * Errors:
 *   ENOMEM  not enough memory (e.g. a filler() call failed).
 *
 * @param path    path to the directory.
 * @param buf     buffer that receives the result.
 * @param filler  function that needs to be called for each directory entry.
 *                Pass 0 as offset (4th argument). 3rd argument can be NULL.
 * @param offset  unused.
 * @param fi      unused.
 * @return        0 on success; -errno on error.
 */
static int a1fs_readdir(const char *path, void *buf, fuse_fill_dir_t filler,
                        off_t offset, struct fuse_file_info *fi)
{
	(void)offset;// unused
	(void)fi;// unused
	//TODO: lookup the directory inode for given path and iterate through its
	// directory entries
	fs_ctx *fs = get_fs();
	struct a1fs_inode *inode;
	find_inode(path, &inode);
	for(unsigned int i = 0; i < 10; i++){
		if(i < inode->extent_count){
			if(filler_entry_in_ext(inode->extents[i], buf, filler)) return -ENOMEM;
		}
	}
	if(inode->extent_count > 10){
		struct a1fs_extent *ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE;
		for(unsigned int j = 0; j < inode->extent_count - 10; j++){
			if(filler_entry_in_ext(ext[j], buf, filler)) return -ENOMEM;
		}
	}
	return 0;
}


/**
 * Create a directory.
 *
 * Implements the mkdir() system call.
 *
 * NOTE: the mode argument may not have the type specification bits set, i.e.
 * S_ISDIR(mode) can be false. To obtain the correct directory type bits use
 * "mode | S_IFDIR".
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" doesn't exist.
 *   The parent directory of "path" exists and is a directory.
 *   "path" and its components are not too long.
 *
 * Errors:
 *   ENOMEM  not enough memory (e.g. a malloc() call failed).
 *   ENOSPC  not enough free space in the file system.
 *
 * @param path  path to the directory to create.
 * @param mode  file mode bits.
 * @return      0 on success; -errno on error.
 */
static int a1fs_mkdir(const char *path, mode_t mode)
{
	fs_ctx *fs = get_fs();

	//TODO: create a directory at given path with given mode
	a1fs_superblock *sb = fs->image;
	if(sb->free_inode_count < 1 || sb->free_block_count < 1) return -ENOSPC;
	char *last_delim = strrchr(path, '/'); //last occurrence
	char *dir_name = last_delim + 1;
	char parent[NAME_MAX];
	strncpy(parent, path, last_delim - path + 1);
	parent[last_delim - path + 1] = '\0';
	struct a1fs_inode *pinode;
	find_inode(parent, &pinode);
	//update parent info
	struct a1fs_dentry *entry;
	if(add_dentry(pinode, &entry)) return -ENOSPC;
	pinode->links += 1;
	pinode->size += sizeof(struct a1fs_dentry);
	//set inode
	int ino_num = set_first_available(1);
	entry->ino = ino_num;
	strcpy(entry->name, dir_name);
	struct a1fs_inode *cinode = fs->image + sb->inode_table * A1FS_BLOCK_SIZE + sizeof(struct a1fs_inode) * ino_num;
	cinode->mode = (mode_t) (mode | S_IFDIR);
	cinode->links = 2;
	cinode->size = 2 * sizeof(struct a1fs_dentry);
	clock_gettime(CLOCK_REALTIME, &(cinode->mtime));
	clock_gettime(CLOCK_REALTIME, &(pinode->mtime));
	cinode->extent_count = 1;
	cinode->extent_blk = 0;
	int blk_num = set_first_available(0);
	cinode->extents[0].start = blk_num;
	cinode->extents[0].count = 1;
	for(int i = 1; i < 10; i++){
		cinode->extents[i].start = 0;
		cinode->extents[i].count = 0;
	}
	//directory entries
	struct a1fs_dentry *first_dentry = fs->image + A1FS_BLOCK_SIZE * blk_num;
	for(unsigned int i = 0; i < A1FS_BLOCK_SIZE / sizeof(struct a1fs_dentry); i++){
		struct a1fs_dentry *dentry = first_dentry + i;
		//indicates unused entry
		dentry->ino = 0; 
		memset(dentry->name, 0, NAME_MAX); 
	}
	return 0;
}


/**
 * Remove a directory.
 *
 * Implements the rmdir() system call.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a directory.
 *
 * Errors:
 *   ENOTEMPTY  the directory is not empty.
 *
 * @param path  path to the directory to remove.
 * @return      0 on success; -errno on error.
 */
static int a1fs_rmdir(const char *path)
{
	fs_ctx *fs = get_fs();

	//TODO: remove the directory at given path (only if it's empty)
	a1fs_superblock *sb = fs->image;
	struct a1fs_inode *cinode;
	find_inode(path, &cinode);
	if(cinode->size > 2 * sizeof(struct a1fs_dentry)) return -ENOTEMPTY;
	char *last_delim = strrchr(path, '/'); //last occurrence
	char *dir_name = last_delim + 1;
	char parent[NAME_MAX];
	strncpy(parent, path, last_delim - path + 1);
	parent[last_delim - path + 1] = '\0';
	struct a1fs_inode *pinode;
	find_inode(parent, &pinode);
	int ino_num = delete_dentry_in_ino(pinode, dir_name);
	//reset block and set block bitmap
	reset_blocks_in_ino(cinode, 0);
	//set inode bitmap
	unsigned char *inode_bitmap = fs->image + sb->inode_bitmap * A1FS_BLOCK_SIZE;
	int byte = ino_num / 8;
	int bit = ino_num % 8;
	inode_bitmap[byte] = inode_bitmap[byte] & ~(1 << bit);
	sb->free_inode_count += 1;
	//update parent inode
	pinode->links -= 1;
	clock_gettime(CLOCK_REALTIME, &(pinode->mtime));

	memset(cinode, 0, sizeof(a1fs_inode));
	return 0;
}


/**
 * Create a file.
 *
 * Implements the open()/creat() system call.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" doesn't exist.
 *   The parent directory of "path" exists and is a directory.
 *   "path" and its components are not too long.
 *
 * Errors:
 *   ENOMEM  not enough memory (e.g. a malloc() call failed).
 *   ENOSPC  not enough free space in the file system.
 *
 * @param path  path to the file to create.
 * @param mode  file mode bits.
 * @param fi    unused.
 * @return      0 on success; -errno on error.
 */
static int a1fs_create(const char *path, mode_t mode, struct fuse_file_info *fi)
{
	(void)fi;// unused
	assert(S_ISREG(mode));
	fs_ctx *fs = get_fs();
	//TODO: create a file at given path with given mode
	a1fs_superblock *sb = fs->image;
	if(sb->free_inode_count < 1) return -ENOSPC; //check free block?
	char *last_delim = strrchr(path, '/'); //last occurrence
	char *dir_name = last_delim + 1;
	char parent[NAME_MAX];
	strncpy(parent, path, last_delim - path + 1);
	parent[last_delim - path + 1] = '\0';
	struct a1fs_inode *pinode;
	find_inode(parent, &pinode);
	//update parent info
	struct a1fs_dentry *entry;
	if(add_dentry(pinode, &entry)) return -ENOSPC;
	pinode->size += sizeof(struct a1fs_dentry);
	//set inode
	int ino_num = set_first_available(1);
	entry->ino = ino_num;
	strcpy(entry->name, dir_name);
	struct a1fs_inode *cinode = fs->image + sb->inode_table * A1FS_BLOCK_SIZE + sizeof(struct a1fs_inode) * ino_num;
	cinode->mode = (mode_t)(mode | S_IFREG);
	cinode->links = 1;
	cinode->size = 0;
	clock_gettime(CLOCK_REALTIME, &(pinode->mtime));
	cinode->extent_count = 0;
	cinode->extent_blk = 0;
	for(int i = 0; i < 10; i++){
		cinode->extents[i].start = 0;
		cinode->extents[i].count = 0;
	}
	return 0;
}

/**
 * Remove a file.
 *
 * Implements the unlink() system call.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a file.
 *
 * @param path  path to the file to remove.
 * @return      0 on success; -errno on error.
 */
static int a1fs_unlink(const char *path)
{
	fs_ctx *fs = get_fs();

	//TODO: remove the file at given path
	a1fs_superblock *sb = fs->image;
	struct a1fs_inode *cinode;
	find_inode(path, &cinode);
	char *last_delim = strrchr(path, '/'); //last occurrence
	char *dir_name = last_delim + 1;
	char parent[NAME_MAX];
	strncpy(parent, path, last_delim - path + 1);
	parent[last_delim - path + 1] = '\0';
	struct a1fs_inode *pinode;
	find_inode(parent, &pinode);
	int ino_num = delete_dentry_in_ino(pinode, dir_name);
	//reset block and set block bitmap
	reset_blocks_in_ino(cinode, 0);
	//set inode bitmap
	unsigned char *inode_bitmap = fs->image + sb->inode_bitmap * A1FS_BLOCK_SIZE;
	int byte = ino_num / 8;
	int bit = ino_num % 8;
	inode_bitmap[byte] = inode_bitmap[byte] & ~(1 << bit);
	sb->free_inode_count += 1;
	//update parent inode
	clock_gettime(CLOCK_REALTIME, &(pinode->mtime));

	memset(cinode, 0, sizeof(a1fs_inode));
	return 0;
}

/**
 * Rename a file or directory.
 *
 * Implements the rename() system call. See "man 2 rename" for details.
 * If the destination file (directory) already exists, it must be replaced with
 * the source file (directory). Existing destination can be replaced if it's a
 * file or an empty directory.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "from" exists.
 *   The parent directory of "to" exists and is a directory.
 *   If "from" is a file and "to" exists, then "to" is also a file.
 *   If "from" is a directory and "to" exists, then "to" is also a directory.
 *
 * Errors:
 *   ENOMEM     not enough memory (e.g. a malloc() call failed).
 *   ENOTEMPTY  destination is a non-empty directory.
 *   ENOSPC     not enough free space in the file system.
 *
 * @param from  original file path.
 * @param to    new file path.
 * @return      0 on success; -errno on error.
 */
static int a1fs_rename(const char *from, const char *to)
{
	fs_ctx *fs = get_fs();
	//TODO: move the inode (file or directory) at given source path to the
	// destination path, according to the description above
	a1fs_superblock *sb = fs->image;
	
	//from
	char *f_last_delim = strrchr(from, '/'); //last occurrence
	char *f_name = f_last_delim + 1;
	char f_parent[NAME_MAX];
	strncpy(f_parent, from, f_last_delim - from + 1);
	f_parent[f_last_delim - from + 1] = '\0';
	struct a1fs_inode *f_pinode;
	find_inode(f_parent, &f_pinode);
	struct a1fs_inode *f_cinode;
	find_inode(from, &f_cinode);
	//to
	char *t_last_delim = strrchr(to, '/'); //last occurrence
	char *t_name = t_last_delim + 1;
	char t_parent[NAME_MAX];
	strncpy(t_parent, to, t_last_delim - to + 1);
	t_parent[t_last_delim - to + 1] = '\0';
	struct a1fs_inode *t_pinode;
	find_inode(t_parent, &t_pinode);
	struct a1fs_inode *t_cinode;
	struct a1fs_dentry *f_dentry;
	struct a1fs_dentry *t_dentry;
	//get from dentry
	f_dentry = find_entry_in_ino(f_pinode, f_name);
	//if to exists
	if(find_inode(to, &t_cinode) == 0){
		if(S_ISDIR(t_cinode->mode) && t_cinode->size > 2 * sizeof(struct a1fs_dentry)) return -ENOTEMPTY;
		reset_blocks_in_ino(t_cinode, 0);
		t_dentry = find_entry_in_ino(t_pinode, t_name);
		unsigned char *inode_bitmap = fs->image + sb->inode_bitmap * A1FS_BLOCK_SIZE;
		memset(t_cinode, 0, sizeof(a1fs_inode));
		int byte = t_dentry->ino / 8;
		int bit = t_dentry->ino % 8;
		inode_bitmap[byte] = inode_bitmap[byte] & ~(1 << bit);
		sb->free_inode_count += 1;
		t_dentry->ino = f_dentry->ino;
	}
	else{
		if(add_dentry(t_pinode, &t_dentry)) return -ENOSPC;
		t_pinode->size += sizeof(struct a1fs_dentry);
		t_dentry->ino = f_dentry->ino;
		strcpy(t_dentry->name, t_name);
		if(S_ISDIR(f_cinode->mode)) t_pinode->links += 1;
	}
	delete_dentry_in_ino(f_pinode, f_name);
	if(S_ISDIR(f_cinode->mode)) f_pinode->links -= 1;
	clock_gettime(CLOCK_REALTIME, &(f_pinode->mtime));
	clock_gettime(CLOCK_REALTIME, &(f_cinode->mtime));
	clock_gettime(CLOCK_REALTIME, &(t_pinode->mtime));
	return 0;
}


/**
 * Change the access and modification times of a file or directory.
 *
 * Implements the utimensat() system call. See "man 2 utimensat" for details.
 *
 * NOTE: You only have to implement the setting of modification time (mtime).
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists.
 *
 * @param path  path to the file or directory.
 * @param tv    timestamps array. See "man 2 utimensat" for details.
 * @return      0 on success; -errno on failure.
 */
static int a1fs_utimens(const char *path, const struct timespec tv[2])
{

	//TODO: update the modification timestamp (mtime) in the inode for given
	// path with either the time passed as argument or the current time,
	// according to the utimensat man page
	
	struct a1fs_inode *inode;
	find_inode(path, &inode);
	if(tv){
		inode->mtime = tv[1];
	}
	else{
		if(clock_gettime(CLOCK_REALTIME, &(inode->mtime)) == -1){
			return -errno;
		}
	}
	return 0;
}

/**
 * Change the size of a file.
 *
 * Implements the truncate() system call. Supports both extending and shrinking.
 * If the file is extended, future reads from the new uninitialized range must
 * return ranges filled with zeros.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a file.
 *
 * Errors:
 *   ENOMEM  not enough memory (e.g. a malloc() call failed).
 *   ENOSPC  not enough free space in the file system.
 *
 * @param path  path to the file to set the size.
 * @param size  new file size in bytes.
 * @return      0 on success; -errno on error.
 */
static int a1fs_truncate(const char *path, off_t size)
{
	//TODO: set new file size, possibly "zeroing out" the uninitialized range
	struct a1fs_inode *inode;
	find_inode(path, &inode);
	if(size > (int)(inode->size)){
		if(inode->extent_count == 0){ //no extents: create one for the file
			int start;
			if((start = set_first_available(0)) == 0) return -ENOSPC;
			inode->extents[0].start = start;
			inode->extents[0].count = 1;
			inode->extent_count = 1;
		}
		if(expand_file(inode, size - (int)(inode->size))) return -ENOSPC;
	}
	else if(size < (int)(inode->size)){
		shrink_file(inode, size);
		inode->size = size;
	}
	return 0;
}


/**
 * Read data from a file.
 *
 * Implements the pread() system call. Should return exactly the number of bytes
 * requested except on EOF (end of file) or error, otherwise the rest of the
 * data will be substituted with zeros. Reads from file ranges that have not
 * been written to must return ranges filled with zeros.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a file.
 *
 * @param path    path to the file to read from.
 * @param buf     pointer to the buffer that receives the data.
 * @param size    buffer size (number of bytes requested).
 * @param offset  offset from the beginning of the file to read from.
 * @param fi      unused.
 * @return        number of bytes read on success; 0 if offset is beyond EOF;
 *                -errno on error.
 */
static int a1fs_read(const char *path, char *buf, size_t size, off_t offset,
                     struct fuse_file_info *fi)
{
	(void)fi;// unused
	fs_ctx *fs = get_fs();

	//TODO: read data from the file at given offset into the buffer
	struct a1fs_inode *inode;
	find_inode(path, &inode);
	if(offset > (int)inode->size){
		return 0;
	}
	unsigned int ext_num = 0;
	int blk_num = 0;
	int remainder = 0;
	if(find_ext_by_offset(inode, offset, &ext_num, &blk_num, &remainder) == 1) return 0;
	unsigned int remain = size;
	struct a1fs_extent *ext;
	int count = 0;
	//get the first extent to read
	if(ext_num <= 9) ext = inode->extents + ext_num;
	else ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE + (ext_num-10) * sizeof(struct a1fs_extent);
	unsigned int size_left = (ext->count - blk_num) * A1FS_BLOCK_SIZE - remainder;
	//read ext
	if(remain >= size_left){
		memcpy(buf, fs->image + (ext->start + blk_num) * A1FS_BLOCK_SIZE + remainder, size_left);
		remain -= size_left;
		ext_num += 1;
		count += size_left;
	}
	else{
		memcpy(buf, fs->image + (ext->start + blk_num) * A1FS_BLOCK_SIZE + remainder, remain);
		count += remain;
		return count;
	}
	//read the extents after ext
	while(remain > 0){
		if(ext_num >= inode->extent_count){ //already the last ext: no more data to read
			memset(buf + count, 0, remain);
			return count;
		}
		//move to next ext
		if(ext_num <= 9) ext = inode->extents + ext_num;
		else ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE + (ext_num-10) * sizeof(struct a1fs_extent);
		if(remain > ext->count * A1FS_BLOCK_SIZE){
			memcpy(buf + count, fs->image + ext->start * A1FS_BLOCK_SIZE, ext->count * A1FS_BLOCK_SIZE);
			remain -= ext->count * A1FS_BLOCK_SIZE;
			ext_num += 1;
			count += ext->count * A1FS_BLOCK_SIZE;
		}
		else{
			memcpy(buf + count, fs->image + ext->start * A1FS_BLOCK_SIZE, remain);
			remain -= remain;
			count += remain;
		}
	}
	return count;
}

/**
 * Write data to a file.
 *
 * Implements the pwrite() system call. Should return exactly the number of
 * bytes requested except on error. If the offset is beyond EOF (end of file),
 * the file must be extended. If the write creates a "hole" of uninitialized
 * data, future reads from the "hole" must return ranges filled with zeros.
 *
 * Assumptions (already verified by FUSE using getattr() calls):
 *   "path" exists and is a file.
 *
 * @param path    path to the file to write to.
 * @param buf     pointer to the buffer containing the data.
 * @param size    buffer size (number of bytes requested).
 * @param offset  offset from the beginning of the file to write to.
 * @param fi      unused.
 * @return        number of bytes written on success; -errno on error.
 */
static int a1fs_write(const char *path, const char *buf, size_t size,
                      off_t offset, struct fuse_file_info *fi)
{	
	(void)fi;// unused
	fs_ctx *fs = get_fs();

	//TODO: write data from the buffer into the file at given offset, possibly
	// "zeroing out" the uninitialized range
	struct a1fs_inode *inode;
	find_inode(path, &inode);
	a1fs_truncate(path, offset + size);

	unsigned int ext_num = 0;
	int blk_num = 0;
	int remainder = 0;
	if(find_ext_by_offset(inode, offset, &ext_num, &blk_num, &remainder) == 1) return 0;
	unsigned int remain = size;
	struct a1fs_extent *ext;
	int count = 0;
	//get the first extent to write on
	if(ext_num <= 9) ext = inode->extents + ext_num;
	else ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE + (ext_num-10) * sizeof(struct a1fs_extent);
	unsigned int size_left = (ext->count - blk_num) * A1FS_BLOCK_SIZE - remainder;
	//write on ext
	if(remain >= size_left){
		memcpy(fs->image + (ext->start + blk_num) * A1FS_BLOCK_SIZE + remainder, buf, size_left);
		remain -= size_left;
		ext_num += 1;
		count += size_left;
	}
	else{
		memcpy(fs->image + (ext->start + blk_num) * A1FS_BLOCK_SIZE + remainder, buf, remain);
		return remain;
	}
	//write on all extents after ext
	while(remain > 0){
		if(ext_num >= inode->extent_count){//already the last extent: ideally wont run into this block since expanded before
			return count;
		}
		if(ext_num <= 9) ext = inode->extents + ext_num;
		else ext = fs->image + inode->extent_blk * A1FS_BLOCK_SIZE + (ext_num-10) * sizeof(struct a1fs_extent);
		if(remain >= ext->count * A1FS_BLOCK_SIZE){
			memcpy(fs->image + ext->start * A1FS_BLOCK_SIZE, buf + count, ext->count * A1FS_BLOCK_SIZE);
			remain -= ext->count * A1FS_BLOCK_SIZE;
			ext_num += 1;
			count += ext->count * A1FS_BLOCK_SIZE;
		}
		else{
			memcpy(fs->image + ext->start * A1FS_BLOCK_SIZE, buf + count, remain);
			remain -= remain;
			count += remain;
		}
	}
	return count;
}


static struct fuse_operations a1fs_ops = {
	.destroy  = a1fs_destroy,
	.statfs   = a1fs_statfs,
	.getattr  = a1fs_getattr,
	.readdir  = a1fs_readdir,
	.mkdir    = a1fs_mkdir,
	.rmdir    = a1fs_rmdir,
	.create   = a1fs_create,
	.unlink   = a1fs_unlink,
	.rename   = a1fs_rename,
	.utimens  = a1fs_utimens,
	.truncate = a1fs_truncate,
	.read     = a1fs_read,
	.write    = a1fs_write,
};

int main(int argc, char *argv[])
{	
	a1fs_opts opts = {0};// defaults are all 0
	struct fuse_args args = FUSE_ARGS_INIT(argc, argv);
	if (!a1fs_opt_parse(&args, &opts)) return 1;

	fs_ctx fs = {0};
	if (!a1fs_init(&fs, &opts)) {
		fprintf(stderr, "Failed to mount the file system\n");
		return 1;
	}

	return fuse_main(args.argc, args.argv, &a1fs_ops, &fs);
}
