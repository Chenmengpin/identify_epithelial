#!/bin/bash

nci_dir="/g/data1a/ku3/jt3341/projects/single_cell/identify_epithelial"
wp_dir="/share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial"
mdss_dir="jt3341/projects/single_cell/identify_epithelial"
ssh jt3341@raijin.nci.org.au "mkdir -p $nci_dir/scripts; mdss mkdir -p $mdss_dir"

rsync -avPS $wp_dir/scripts/* jt3341@raijin.nci.org.au:$nci_dir/scripts

ssh jt3341@raijin.nci.org.au "cd $nci_dir; tar -zcvf scripts.tar.gz scripts; \
	mdss put $nci_dir/scripts.tar.gz $mdss_dir; rm -r $nci_dir/scripts.tar.gz $nci_dir/scripts"

cd $wp_dir/scripts
git add -A
git commit -m "regular"
git push