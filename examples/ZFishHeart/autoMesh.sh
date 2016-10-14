#!/bin/bash
FILES=FolderToMyData/*
for f in $FILES
do
  echo "Processing $f file..."
  # take action on each file. $f store current file name
  
  meshlabserver -i $f -o ${f%.obj}_mesh.ply -s heartFinal1.mlx -om vn

  # ${f%.obj}_mesh.obj
done
