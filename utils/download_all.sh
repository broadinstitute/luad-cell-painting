#!/bin/bash

# Download data from https://registry.opendata.aws/cell-painting-image-collection/
# WARNING!. This script will download metadata, raw tif images and SQLite
# features. It will use ~600 GB.
read -p "You will download and override your data(~600GB). Continue (yes/no)?" CONT
echo $CONT
if [ "$CONT" != "yes" ]; then
  echo "No changes applied.";
  exit 0
fi
rootpath=$(dirname $0)/..
# Metadata is now included in the repo.
# echo "Downloading Metadata data..."
# aws s3 sync s3://cytodata/datasets/LUAD-BBBC043-Caicedo/metadata/LUAD-BBBC043-Caicedo/ $rootpath/inputs/metadata/cytodata --no-sign-request
echo "Downloading TIF Images data..."
aws s3 sync s3://cytodata/datasets/LUAD-BBBC043-Caicedo/images/LUAD-BBBC043-Caicedo/ $rootpath/inputs/images/ --no-sign-request
echo "Downloading SQLite data..."
plates=(52649 52651 52652 52653 52654 52657 52662 52663 52664 52665 52666 52671 52672 52673 52674 52675)
mkdir -p inputs/locations/sqlite
for plate in ${plates[@]}; do
    aws s3 cp s3://cytodata/datasets/LUAD-BBBC043-Caicedo/profiles_cp/LUAD-BBBC043-Caicedo/$plate/$plate.sqlite inputs/locations/sqlite/ --no-sign-request
done
