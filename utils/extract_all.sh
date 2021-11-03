#!/bin/bash
plates=(52649 52651 52652 52653 52654 52657 52662 52663 52664 52665 52666 52671 52672 52673 52674 52675)
prefixpath=$(dirname $0)
for plate in ${plates[@]}; do
    python $prefixpath/extract_locations.py inputs/locations/sqlite/$plate.sqlite inputs/locations/$plate/
done
