# Construction of mitochondrial reference genomes of the Presbytis femoralis group using multiple sequence alignment for new species identification
## ZB4171-Mandolin: Ng Jing Ting, Quek Zhi Heng, and Tay Hui Yi

This project was done under ZB4171, a module from the National University of Singapore. 

Our project explores prior work done by Meier et. al (2020). Contrary to what was done in their work, we suggest the use of all available samples for P. f. percura and P. f. presbytis to construct multiple proposed reference genomes for each sample before using multiple sequence alignment to produce a final reference genome that chooses the most conserved sample instead of relying solely on coverage. A comparison between the results of this method and the paper could shed light on potential sample biases and serve to validate the results of the study.

This repository contains a the dataset and code files required to replicate the bioinformatics analysis done by Meier et. al. This workflow is conducted on AWS, using Nextflow. The detailed set-up may be found below. 

## EC2 specificiations and set up
- Recommended instance type: c5.4xlarge
- Recommended storage size: >= 250 GB
  
In the case where an extra volume is required, you may need to create and mount the new volume: https://devopscube.com/mount-ebs-volume-ec2-instance/

## Downloading conda
- Refer to the set-up here: https://docs.conda.io/projects/miniconda/en/latest/
- To restart terminal session without quitting, type “source ~/.bashrc

## Installing Nextflow
- Refer to the set-up here: https://anaconda.org/bioconda/nextflow
- To check if Nextflow is downloaded successfully, run $nextflow run hello.

## Installing Docker
- Refer to the set-up here: https://www.cyberciti.biz/faq/how-to-install-docker-on-amazon-linux-2/
- Check if docker can pull image : $docker pull hello-world
  
Potential problems encountered: Configure docker (if permission denied in above)
- sudo usermod -a -G docker $USER
- reboot

## Downloading fastq files using SRA tools
- Refer to the set-up here: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
- Follow steps 1, 2, 3, and 8 (minimally).
  
To bypass exporting the path every time you enter the machine, do the following 
- vi ~/.bashrc
- Go to the bottom of the file 
- Press “o” to insert a new line
- Add the following command at the bottom of the file 
- export PATH=$PATH:$PWD/sratoolkit.3.0.7-ubuntu64/bin
- Press “esc” and type “:x” to save file
- source ~/.bashrc

Potential problems encountered: aws pkcs7 error when using prefetch
- go to aws console under actions > instance settings > modify instance metadata options > change IMDSv2 from required to optional

Carry out prefetch and fastq-dump (eg):
- $prefetch SRR10277272
- $fastq-dump SRR10277272 --split-3 --skip-technical --gzip

## To create and attach shared volume to your own EC2 instance (note to do this before proceeding to add other volumes):
1. Go to AWS Console and under EC2, click “Volumes” at the side and then “Create Volume”
2. Select “io2” for Volume Type, enable multi-attach and ensure that the availability zone is correct (same as your instance)
3. Click on the link of the volume you just created and copy paste the volume ID of the newly created volume and paste it below once you are done with all the steps
4. Under actions select “Attach Volume”
5. Select your instance that you are running (if you don't see it likely the ec2 instance is not in ap-southeast-1a,change it)
6. For device name, change where appropriate (eg. the last letter of device name) 
7. Once the attaching and all is done, ssh into your ec2 instance
8. Type the command `lsblk` and check that there is a new disk is there (likely it's called nvme1n1)
9. Type the command: ```sudo mkfs -t xfs /dev/nvme1n1```
10. Type the command `sudo file -s /dev/nvme1n1` and check that the output should looks something like:  “SGI XFS filesystem data (blksz 4096, inosz 512, v2 dirs)”
11. Type the command: ```sudo mkdir /shared_drive_{your initials}``` eg.(`sudo mkdir /shared_drive_zh`)
12. Type the command: ```sudo chmod 777 /shared_drive_{your initials}```
13. Type the command: ```sudo mount /dev/nvme1n1 /shared_drive_{your initials}```
14. Type `sudo cp /etc/fstab /etc/fstab.orig` to back up the original file in case there are errors
15. Type `sudo blkid` and copy the UUID of the /dev/nvme1n1 device 
16. Type `sudo vim /etc/fstab`
17. Add the following line to the end of the file: ```UUID={copied UUID from step 11}  /shared_drive_{your initials}  xfs  defaults,nofail  0  2```
18. Save the changes on the file and exit Vim
19. Type the command: `sudo umount /shared_drive_{your initials}`
20. Type the command `sudo mount -a`
21. Reboot your ec2 instance and ssh in again
22. Type `df -h`` and see if the /shared_drive_{your initials} path appears in the last column

## To add other members volumes to your instance:
- You would require the links of the volume IDs of other members
- Repeat steps 1, 2, 3, 9. For each new volume that you attach, the nvme path will change so its likely to be nvme1n2/nvme1n3, for that just change the path accordingly in the respective steps.
- Replace the initials (names) where appropriate

## MITOS Annotation Results

| Accession Number | Sample ID | Species | MITOS Annotation Link |
|:----------------:|:---------:|:-------:|:---------------:|
| SRR10277267 | BLM5 | *Presbytis femoralis femoralis* | [BLM5](http://mitos.bioinf.uni-leipzig.de/result.py?hash=xa4K05aI) |
| SRR10277268 | BLM4 | *Presbytis femoralis femoralis* | [BLM4](http://mitos.bioinf.uni-leipzig.de/result.py?hash=PT2Aj1sw) |
| SRR10277269 | BLM3 | *Presbytis femoralis femoralis* | [BLM3](http://mitos.bioinf.uni-leipzig.de/result.py?hash=Obpy1si2) |
| SRR10277270 | BLM2 | *Presbytis femoralis femoralis* | [BLM2](http://mitos.bioinf.uni-leipzig.de/result.py?hash=d6tbfPru) |
| SRR10277271 | BLM1 | *Presbytis femoralis femoralis* | [BLM1](http://mitos.bioinf.uni-leipzig.de/result.py?hash=PTWwgweN) |
| SRR10277272 | Pres2 | *Presbytis siamensis cf. cana* | [Pres2](http://mitos.bioinf.uni-leipzig.de/result.py?hash=UUpNqQ8k) |
| SRR10277273 | ESBL_8b | *Presbytis femoralis percura* | No reconstructed reference |
| SRR10277274 | ESBL_6a | *Presbytis femoralis percura* | [ESBL_6a](http://mitos.bioinf.uni-leipzig.de/result.py?hash=iOpxEfIY) |
| SRR10277275 | BLM6 | *Presbytis femoralis femoralis* | [BLM6](http://mitos.bioinf.uni-leipzig.de/result.py?hash=tYXWqH1k) |
| SRR10277276 | ESBL_5 | *Presbytis femoralis percura* | [ESBL_5](http://mitos.bioinf.uni-leipzig.de/result.py?hash=HwSBUtGg) |
| SRR10277277 | ESBL_1a | *Presbytis femoralis percura* | [ESBL_1a](http://mitos.bioinf.uni-leipzig.de/result.py?hash=m3il2tFj) |
