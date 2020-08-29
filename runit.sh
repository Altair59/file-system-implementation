echo "====creat and write img===="
truncate -s 65536 img
./mkfs.a1fs -i 16 img
echo "====create folder and mount fs===="
mkdir ~/369
./a1fs img ~/369
echo "====create directory and write something to a file===="
cd ~/369
mkdir d
cd d
touch f
ls
echo csc369 > f
echo assignment 1 >> f
cat f
echo "====unmount and remount===="
cd ~/group_0299/A1b
fusermount -u ~/369
./a1fs img ~/369
echo "====display the contents still exist===="
cd ~/369
ls
cd d
ls
cat f
