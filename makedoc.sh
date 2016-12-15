#!/bin/bash
set -e
echo
echo === Running Javadoc ===
echo
javadoc -sourcepath ./src/ -d ./doc -subpackages . -cp $COATJAVA/lib/clas/coat-libs-3.0-SNAPSHOT.jar -public
echo
echo === Deleting old jar file ===
echo
rm SVTGeometryDoc.jar
echo
echo === Creating new jar file ===
echo
jar vcf SVTGeometryDoc.jar doc/*
