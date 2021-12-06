# nngd

Build nearest-neighbors graph in a ~efficient way.

Compute graph distances + ECTD distances.

Use file knncpp.h from https://github.com/Rookfighter/knn-cpp :
wget https://raw.githubusercontent.com/Rookfighter/knn-cpp/master/include/knncpp.h

## Usage

roxygenize(".") <br>
R CMD INSTALL . <br>
library(nngd) <br>
?nng <br>
?ectd <br>
