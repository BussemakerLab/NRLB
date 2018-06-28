#/!/bin/bash

mkdir temp_bin
javac -cp src/selex.jar -d temp_bin/ src/*/*.java
cp src/selex.jar temp_bin/NRLB.jar
cd temp_bin
jar uf NRLB.jar Jama/*
jar uf NRLB.jar base/*
jar uf NRLB.jar dynamicprogramming/*
jar uf NRLB.jar minimizers/*
jar uf NRLB.jar model/*
jar uf NRLB.jar utils/*

cd ../
mv temp_bin/NRLB.jar bin/
rm -rf temp_bin/
