# numc

## Here's what I did in project 4:

I started by creating the most basic (and most likely inefficient) functions for each operation. This allowed me to fully grasp not only the correctness of the functions themselves, but also how the debugger worked. Many times, I had to use GDB to step through my functions and ouptut the value of a specific variable at every step, and making the functions simplistic helped me through this process.Then I moved on to making these functions more efficient. This involved tasks such as getting rid of function calling and replacing the matrix data mutation by offsets (so instead of set(matrix, i, j), I did matrix(i\*n + j) = x). This made my process of optimizing matrices in step 4 much easier.

I cannot stress enough the importance of gdb for my debugging my functions in step 1. I initially had lots of segfaults at random areas of my functions, but gdb helped me expedite the process of identifying and fixing these segfaults much easier. I also found a trick to help me debug fast. Because I didn't know how to setup valgrind for this project, what i did was use fprintf to debug my variables at each stage of the function. For example, before a segfault diagnosed by gdb, I would place an fprintf on each variable that I supected caused the segfault. This helped me identify where there occured a potential memory pointer error much quicker.

Then, I went on to numc.c. I found that perfecting one function in numc like Matrix_abs and using that result for the other similar ones like Matrix_add, Matrix_sub, Matrix_neg, made my process much easier. This is because I realized the similarity between the four functions. Often times, I only had to change 1-2 lines of code when copying one code from one function to the other.

Pow and mul were a different story. I was super confused about how to turn the column-major formatted functions given in lectures to row-major. I initially tried random combinations of code, and wanted to see if any of them worked, but none of them did. Instead, I pulled out a piece of paper, drew 2 10x10 matrices, and walked through each loop of the provided solution in lecture to see what they actually did. Then, with the patterns I noticed, I tried to reproduce it within my row-major formatted matrices and it ended up working.

For simd, multi-threading, and unrolling, this was a matter of labs: I basically went through my answers for each lab and tried to reproduce the results. The hardest part for me was the tail cases, but in the end, I managed to work around it my drawing similarities and differences between the lab code and my own.

One area that I should mention caused alot of confusion for me was the errors that we were told to output when, for example, the user inputed a matrix with the wrong dimension or when the user inputed value that was not a double. I had trouble understanding how to use the inbuild util functions such as PyObject_TypeCheck, PyList_Check etc. But by far the most useful tools for me to learn and understand how to find these errors were the prewritten functions such as Matrix61c_init.

