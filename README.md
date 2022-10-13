Building: simply do

$ g++ -o bitcube main.cpp

in the unzipped directory.

Usage: since this executable is for encrypting/decrypting the signed decimal numbers, for encrypting the QDCT coefficients within JPEG compression, it must be exploited along with the modified libjpeg-turbo as follows:

$ cjpeg -wqdct tmp1.y tmp1.u tmp1.v -sample 2x2 -q 70 -outfile out.jpg in.ppm
$ bitcube tmp1.y tmp1.u tmp1.v tmp2.y tmp2.u tmp2.v
$ cjpeg -rqdct tmp2.y tmp2.u tmp2.v -sample 2x2 -q 70 -outfile out.jpg in.ppm

For the decoder side, similarly do

$ djpeg -wqdct tmp1.y tmp1.u tmp1.v -outfile out.ppm in.jpg
$ bitcube tmp1.y tmp1.u tmp1.v tmp2.y tmp2.u tmp2.v
$ djpeg -rqdct tmp2.y tmp2.u tmp2.v -outfile out.ppm in.jpg


