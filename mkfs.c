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
 * CSC369 Assignment 1 - a1fs formatting tool.
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include <time.h>

#include "a1fs.h"
#include "map.h"


/** Command line options. */
typedef struct mkfs_opts {
	/** File system image file path. */
	const char *img_path;
	/** Number of inodes. */
	size_t n_inodes;

	/** Print help and exit. */
	bool help;
	/** Overwrite existing file system. */
	bool force;
	/** Sync memory-mapped image file contents to disk. */
	bool sync;
	/** Verbose output. If false, the program must only print errors. */
	bool verbose;
	/** Zero out image contents. */
	bool zero;

} mkfs_opts;

static const char *help_str = "\
Usage: %s options image\n\
\n\
Format the image file into a1fs file system. The file must exist and\n\
its size must be a multiple of a1fs block size - %zu bytes.\n\
\n\
Options:\n\
    -i num  number of inodes; required argument\n\
    -h      print help and exit\n\
    -f      force format - overwrite existing a1fs file system\n\
    -s      sync image file contents to disk\n\
    -v      verbose output\n\
    -z      zero out image contents\n\
";

static void print_help(FILE *f, const char *progname)
{
	fprintf(f, help_str, progname, A1FS_BLOCK_SIZE);
}


static bool parse_args(int argc, char *argv[], mkfs_opts *opts)
{
	char o;
	while ((o = getopt(argc, argv, "i:hfsvz")) != -1) {
		switch (o) {
			case 'i': opts->n_inodes = strtoul(optarg, NULL, 10); break;

			case 'h': opts->help    = true; return true;// skip other arguments
			case 'f': opts->force   = true; break;
			case 's': opts->sync    = true; break;
			case 'v': opts->verbose = true; break;
			case 'z': opts->zero    = true; break;

			case '?': return false;
			default : assert(false);
		}
	}

	if (optind >= argc) {
		fprintf(stderr, "Missing image path\n");
		return false;
	}
	opts->img_path = argv[optind];

	if (opts->n_inodes == 0) {
		fprintf(stderr, "Missing or invalid number of inodes\n");
		return false;
	}
	return true;
}


/** Determine if the image has already been formatted into a1fs. */
static bool a1fs_is_present(void *image)
{
	//TODO: check if the image already contains a valid a1fs superblock
	struct a1fs_superblock *sb = image;
	if(sb->magic == A1FS_MAGIC){
		return true;
	}
	return false;
}


/**
 * Set the bitmap bit to 1 according to num.
 * 
 * @param image  pointer to the start of the image.
 * @param ibmp   1 if the bitmap to be set is inode bitmap, 0 if the bitmap to be set is block bitmap.
 * @param num    the position of the bit to be set, starting from 0.
 */
void set_bitmap_to_1(void *image, int ibmp, int num){
	struct a1fs_superblock *sb = image;
	unsigned char *start;
	if(ibmp) start = image + sb->inode_bitmap * A1FS_BLOCK_SIZE;
	else start = image + sb->block_bitmap * A1FS_BLOCK_SIZE;	
	int byte = num / 8;
	int bit = num % 8;
	start[byte] = start[byte] | (1 << bit);
}


/**
 * Format the image into a1fs.
 *
 * NOTE: Must update mtime of the root directory.
 *
 * @param image  pointer to the start of the image.
 * @param size   image size in bytes.
 * @param opts   command line options.
 * @return       true on success;
 *               false on error, e.g. options are invalid for given image size.
 */
static bool mkfs(void *image, size_t size, mkfs_opts *opts)
{
	//TODO: initialize the superblock and create an empty root directory
	memset(image, 0, size);
	int num_inodes = opts->n_inodes;
	int num_blocks = size / A1FS_BLOCK_SIZE;
	struct a1fs_superblock *sb = image;
	sb->magic = (uint64_t)A1FS_MAGIC;
	sb->size = (uint64_t)size;
	sb->inode_count = (uint32_t)(num_inodes);
	sb->free_inode_count = (uint32_t)(num_inodes - 1);
	sb->inode_size = (uint32_t)sizeof(struct a1fs_inode);
	sb->block_count = (uint32_t)num_blocks;
	//Inode bitmap
	//unsigned char *inode_bitmap = image + A1FS_BLOCK_SIZE;
	int num_ibmp_block = num_inodes / (A1FS_BLOCK_SIZE * 8);
	if(num_inodes % (A1FS_BLOCK_SIZE * 8) > 0) num_ibmp_block += 1;
	sb->inode_bitmap = (a1fs_blk_t)1;
	sb->block_bitmap = (a1fs_blk_t)(sb->inode_bitmap + num_ibmp_block);
	set_bitmap_to_1(image, 1, 0);

	//Counting number of inode table blocks.
	int num_table_blk = num_inodes * sb->inode_size / A1FS_BLOCK_SIZE;
	if((num_inodes * sb->inode_size) % A1FS_BLOCK_SIZE > 0) num_table_blk += 1;
	//Block bitmap
	//unsigned char *block_bitmap = inode_bitmap + A1FS_BLOCK_SIZE * num_ibmp_block;
	int num_bbmp_block = num_blocks / (A1FS_BLOCK_SIZE * 8);
	if(num_blocks % (A1FS_BLOCK_SIZE * 8) > 0) num_bbmp_block += 1;
	int num_taken = 1 + num_ibmp_block + num_table_blk + num_bbmp_block + 1;
	if(num_taken > num_blocks) return false;
	for(int bit = 0; bit < num_taken; bit++){
		set_bitmap_to_1(image, 0, bit);
	}
	
	sb->inode_table = (a1fs_blk_t)(sb->block_bitmap + num_bbmp_block);
	sb->free_block_count = (uint32_t)(num_blocks - num_taken);
	sb->data_start = sb->inode_table + num_table_blk;
	//Inode table and inode for root
	//struct a1fs_inode *inode_table = block_bitmap + A1FS_BLOCK_SIZE * num_bbmp_block;
	struct a1fs_inode *first_inode = image + A1FS_BLOCK_SIZE * sb->inode_table;
	first_inode->mode = (mode_t)(S_IFDIR | 0777);
	first_inode->links = (uint32_t)2;
	first_inode->size = (uint64_t)(2 * sizeof(a1fs_dentry));
	clock_gettime(CLOCK_REALTIME, &(first_inode->mtime));
	first_inode->extent_count = (uint32_t)1;
	first_inode->extent_blk = (uint32_t)0;
	first_inode->extents[0].start = (a1fs_blk_t)(sb->inode_table + num_table_blk);
	first_inode->extents[0].count = (a1fs_blk_t)1;
	for(int i = 1; i < 10; i++){
		first_inode->extents[i].start = 0;
		first_inode->extents[i].count = 0;
	}
	//directory entries of root
	struct a1fs_dentry *first_dentry = image + A1FS_BLOCK_SIZE * sb->data_start;
	for(unsigned int i = 0; i < A1FS_BLOCK_SIZE / sizeof(struct a1fs_dentry); i++){
		struct a1fs_dentry *dentry = first_dentry + i;
		//indicates unused entry
		dentry->ino = 0;
		memset(dentry->name, 0, NAME_MAX); 
	}
	return true;
}


int main(int argc, char *argv[])
{
	mkfs_opts opts = {0};// defaults are all 0
	if (!parse_args(argc, argv, &opts)) {
		// Invalid arguments, print help to stderr
		print_help(stderr, argv[0]);
		return 1;
	}
	if (opts.help) {
		// Help requested, print it to stdout
		print_help(stdout, argv[0]);
		return 0;
	}

	// Map image file into memory
	size_t size;
	void *image = map_file(opts.img_path, A1FS_BLOCK_SIZE, &size);
	if (image == NULL) return 1;

	// Check if overwriting existing file system
	int ret = 1;
	if (!opts.force && a1fs_is_present(image)) {
		fprintf(stderr, "Image already contains a1fs; use -f to overwrite\n");
		goto end;
	}

	if (opts.zero) memset(image, 0, size);
	if (!mkfs(image, size, &opts)) {
		fprintf(stderr, "Failed to format the image\n");
		goto end;
	}

	// Sync to disk if requested
	if (opts.sync && (msync(image, size, MS_SYNC) < 0)) {
		perror("msync");
		goto end;
	}

	if(opts.verbose){
		printf("An implementation of extent-based file system.\nChanges made from the descriptions in the proposal:\n1. Removed ctime in inode.\n2. Added one more extent in the extent array of each inode.\n3. Added padding in the inode.\n");
	}

	ret = 0;
end:
	munmap(image, size);
	return ret;
}