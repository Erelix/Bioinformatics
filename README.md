# Bioinformatics

py -3.12 lab1.py


KODONAS

(M1:0.2662,(B3:0.1994,(B1:0.1687,(B4:0.3729,M3:0.2990):0.0401):0.0937):0.0479,(B2:0.1433,(M2:0.1107,M4:0.3044):0.1131):0.1704);

Weight matix is of the form W=1/D^0.000000 

Tree reconstruction method - Weighted least-squares method MW (global optimization)


TREE METRIC (ADDITIVE DISTANCE) MATRIX (AD)

B1                            0.000000 0.062391 0.046181 0.058174 0.057651 0.070446 0.050778 0.089810 
B2                            0.062391 0.000000 0.056091 0.086822 0.057988 0.036709 0.079426 0.056073 
B3                            0.046181 0.056091 0.000000 0.070613 0.051352 0.064147 0.063217 0.083511 
B4                            0.058174 0.086822 0.070613 0.000000 0.082083 0.094878 0.067192 0.114242 
M1                            0.057651 0.057988 0.051352 0.082083 0.000000 0.066043 0.074687 0.085407 
M2                            0.070446 0.036709 0.064147 0.094878 0.066043 0.000000 0.087482 0.041510 
M3                            0.050778 0.079426 0.063217 0.067192 0.074687 0.087482 0.000000 0.106846 
M4                            0.089810 0.056073 0.083511 0.114242 0.085407 0.041510 0.106846 0.000000 


THE FOLLOWING STATISTICS ARE AVAILABLE FOR 
A GIVEN DISSIMILARITY (D) AND AN OBTAINED TREE METRIC (AD)


Least-squares coefficient  Sum (Dij-ADij)^2  = 0.0003735881
                            i<j 
Average absolute difference  Sum |Dij-ADij|/(n(n-1)/2)  = 0.0028781161
                             i<j 
Maximum absolute difference   Max|Dij-ADij| = 0.0086831878
                              i,j 
Total length of the tree   L  = 0.2329761219


              TREE EDGES WITH THEIR LENGTHS

               9--4               0.037294
              10--9               0.004009
               9--7               0.029898
              11--3               0.019941
              12--5               0.026624
              13--2               0.014327
              14--6               0.011073
               1--10              0.016871
              10--11              0.009369
              11--12              0.004787
              12--13              0.017037
              13--14              0.011309
              14--8               0.030437




DIKODON

(M1:0.9825,(B3:0.8475,(B1:0.8395,(B4:1.3523,M3:1.0172):0.0858):0.1091):0.0638,(B2:0.8005,(M2:0.6709,M4:0.9792):0.1194):0.2010);


Weight matix is of the form W=1/D^0.000000 

Tree reconstruction method - Weighted least-squares method MW (global optimization)


TREE METRIC (ADDITIVE DISTANCE) MATRIX (AD)

B1                            0.000000 0.020138 0.017960 0.022775 0.019948 0.020036 0.019424 0.023119 
B2                            0.020138 0.000000 0.019128 0.026124 0.019840 0.015908 0.022773 0.018991 
B3                            0.017960 0.019128 0.000000 0.023946 0.018937 0.019026 0.020595 0.022109 
B4                            0.022775 0.026124 0.023946 0.000000 0.025933 0.026022 0.023695 0.029105 
M1                            0.019948 0.019840 0.018937 0.025933 0.000000 0.019738 0.022583 0.022821 
M2                            0.020036 0.015908 0.019026 0.026022 0.019738 0.000000 0.022671 0.016502 
M3                            0.019424 0.022773 0.020595 0.023695 0.022583 0.022671 0.000000 0.025754 
M4                            0.023119 0.018991 0.022109 0.029105 0.022821 0.016502 0.025754 0.000000 


THE FOLLOWING STATISTICS ARE AVAILABLE FOR 
A GIVEN DISSIMILARITY (D) AND AN OBTAINED TREE METRIC (AD)


Least-squares coefficient  Sum (Dij-ADij)^2  = 0.0000083402
                            i<j 
Average absolute difference  Sum |Dij-ADij|/(n(n-1)/2)  = 0.0004192914
                             i<j 
Maximum absolute difference   Max|Dij-ADij| = 0.0012051365
                              i,j 
Total length of the tree   L  = 0.0806854125


              TREE EDGES WITH THEIR LENGTHS

               9--7               0.010172
              10--9               0.000858
               9--4               0.013523
              11--3               0.008475
              12--5               0.009825
              13--2               0.008005
              14--6               0.006709
               1--10              0.008395
              10--11              0.001091
              11--12              0.000638
              12--13              0.002010
              13--14              0.001194
              14--8               0.009792