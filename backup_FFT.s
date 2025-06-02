#define STDOUT 0xd0580000

.section .text
.global _start
_start:
    la a0, matrix
    la t0, size
    lw a1, 0(t0)         # Load matrix size into a1
    call bit_reverse

    la a0, matrix
    la t0, size
    lw a1, 0(t0)
    call fft_all_stages

    la a0, matrix
    la t0, size
    lw a1, 0(t0)
    call printToLog
    j _finish

############################################################
# Bit reversal of array based on index
bit_reverse:
    addi sp, sp, -20
    sw ra, 16(sp)
    sw s0, 12(sp)
    sw s1, 8(sp)
    sw s2, 4(sp)
    sw s3, 0(sp)

    mv s0, a0      # base address of matrix
    mv s1, a1      # N

    li t0, 0       # i = 0
loop_i:
    bge t0, s1, done_bit_reverse
    mv a2, t0
    mv a3, s1
    call reverse_bits
    mv t1, a0      # reversed index

    blt t0, t1, swap_elements

next_i:
    addi t0, t0, 1
    j loop_i

swap_elements:
    slli t2, t0, 2
    add t2, s0, t2
    flw ft0, 0(t2)

    slli t3, t1, 2
    add t3, s0, t3
    flw ft1, 0(t3)

    fsw ft1, 0(t2)
    fsw ft0, 0(t3)
    j next_i

done_bit_reverse:
    lw ra, 16(sp)
    lw s0, 12(sp)
    lw s1, 8(sp)
    lw s2, 4(sp)
    lw s3, 0(sp)
    addi sp, sp, 20
    ret

############################################################
# Reverses bits of a2 given log2(a3)
reverse_bits:
    li t0, 0           # result
    mv t1, a2          # input index
    li t2, 0           # bit count

    mv t3, a3          # N
calc_log2_loop:
    srli t3, t3, 1
    beqz t3, done_calc_log2
    addi t2, t2, 1
    j calc_log2_loop
done_calc_log2:

reverse_loop:
    beqz t2, done_reverse
    andi t4, t1, 1
    slli t0, t0, 1
    or t0, t0, t4
    srli t1, t1, 1
    addi t2, t2, -1
    j reverse_loop

done_reverse:
    mv a0, t0
    ret
############################################################
# Perform all FFT stages
# a0 = base address of matrix
# a1 = N (size)
fft_all_stages:
    addi sp, sp, -16
    sw ra, 12(sp)
    sw s0, 8(sp)
    sw s1, 4(sp)
    sw s2, 0(sp)

    mv s0, a0      # base address
    mv s1, a1      # N

    call log2      # a1 = N → a0 = log2(N)
    mv t2, a0      # t2 = total stages

    li t0, 1       # current stage = 1
stage_loop:
    bgt t0, t2, end_stages

    mv a0, s0      # base address
    mv a1, s1      # size N
    mv a2, t0      # current stage
    call fft_stageN

    addi t0, t0, 1
    j stage_loop

end_stages:
    lw ra, 12(sp)
    lw s0, 8(sp)
    lw s1, 4(sp)
    lw s2, 0(sp)
    addi sp, sp, 16
    ret

############################################################
# Performs butterfly operation for a given stage
# a0 = base address, a1 = N, a2 = stage
fft_stageN:
    addi sp, sp, -16
    sw ra, 12(sp)
    sw s0, 8(sp)
    sw s1, 4(sp)
    sw s2, 0(sp)

    mv s0, a0      # base
    mv s1, a1      # N
    mv s2, a2      # stage

    li t0, 1
    sll t1, t0, s2     # m = 2^stage
    srli t2, t1, 1      # half_m = m / 2

    li t3, 0            # k = 0
outer_k_loop:
    bge t3, s1, end_fft_stageN

    li t4, 0            # j = 0
inner_j_loop:
    bge t4, t2, next_k_group

    add t5, t3, t4          # i = k + j
    sll t6, t5, 2
    add a3, s0, t6
    flw ft0, 0(a3)          # x0 = A[i]

    add a4, t5, t2          # i + half_m
    sll a5, a4, 2
    add a6, s0, a5
    flw ft1, 0(a6)         # x1 = A[i + half_m]

    fadd.s ft2, ft0, ft1
    fsub.s ft3, ft0, ft1

    fsw ft2, 0(a3)
    fsw ft3, 0(a6)

    addi t4, t4, 1
    j inner_j_loop

next_k_group:
    add t3, t3, t1          # k += m
    j outer_k_loop

end_fft_stageN:
    lw ra, 12(sp)
    lw s0, 8(sp)
    lw s1, 4(sp)
    lw s2, 0(sp)
    addi sp, sp, 16
    ret


transpose:
    # Prologue
    addi sp, sp, -16
    sw ra, 12(sp)
    sw s0, 8(sp)
    sw s1, 4(sp)
    sw s2, 0(sp)

    mv s0, a0      # s0 = base address of matrix
    mv s1, a1      # s1 = size N

    li t0, 0       # i = 0 (outer loop counter)

outer_loop:
    bge t0, s1, end_transpose  # if i >= N, exit

    addi t1, t0, 1  # j = i + 1 (inner loop counter)

inner_loop:
    bge t1, s1, next_row  # if j >= N, go to next row

    # Compute offsets: (i * N + j) and (j * N + i)
    mul t2, t0, s1  # t2 = i * N
    add t2, t2, t1  # t2 = i * N + j
    slli t2, t2, 2  # t2 *= 4 (float size)

    mul t3, t1, s1  # t3 = j * N
    add t3, t3, t0  # t3 = j * N + i
    slli t3, t3, 2  # t3 *= 4

    add t2, t2, s0  # Address of matrix[i][j]
    add t3, t3, s0  # Address of matrix[j][i]

    # Swap elements
    flw ft0, 0(t2)  # Load matrix[i][j] into ft0
    flw ft1, 0(t3)  # Load matrix[j][i] into ft1
    fsw ft1, 0(t2)  # Store matrix[j][i] into matrix[i][j]
    fsw ft0, 0(t3)  # Store matrix[i][j] into matrix[j][i]

    addi t1, t1, 1  # j++
    j inner_loop

next_row:
    addi t0, t0, 1  # i++
    j outer_loop

end_transpose:
    # Epilogue
    lw ra, 12(sp)
    lw s0, 8(sp)
    lw s1, 4(sp)
    lw s2, 0(sp)
    addi sp, sp, 16
    ret
## END YOU CODE HERE

# Function to print a matrix for debugging purposes
# This function iterates over all elements of a matrix stored in memory.
# Instead of calculating the end address in each iteration, it precomputes 
# the end address (baseAddress + size^2 * 4) to optimize the loop.
# Input:
#   a0: Base address of the matrix
#   a1: Size of matrix
# Clobbers:
#   t0, t1, ft0
printToLog:
    li t0, 0x123                #  Identifiers used for python script to read logs
    li t0, 0x456
    mv a1, a1                   # moving size to get it from log 
    mv t0, a0                   # Copy the base address of the matrix to t0 to avoid modifying a0
    mul t1, a1, a1              # size^2 
    slli  t1, t1, 2             # size^2 * 4 (total size of the matrix in bytes)
    add t1, a0, t1              # Calculate the end address (base address + total size)

    printMatrixLoop:
        bge t0, t1, printMatrixLoopEnd 
        flw ft0, 0(t0)          # Load from array
        addi t0, t0, 4          # increment address by elem size
        j printMatrixLoop
    printMatrixLoopEnd:

    li t0, 0x123                #  Identifiers used for python script to read logs
    li t0, 0x456

    jr ra


# Function: _finish
# VeeR Related function which writes to to_host which stops the simulator
_finish:
    li x3, 0xd0580000
    addi x5, x0, 0xff
    sb x5, 0(x3)
    beq x0, x0, _finish

    .rept 100
        nop
    .endr


.data
## ALL DATA IS DEFINED HERE LIKE MATRIX, CONSTANTS ETC

## DATA DEFINE START
.equ MatrixSize, 8
matrix:
    .float 898.25, -510.5, 49.0, -46.5, -156.75, 448.0, 811.75, -554.25
    .float -393.0, -869.5, -114.0, -0.25, -829.0, 447.25, -321.0, -978.75
    .float 985.25, -377.0, 979.75, -31.0, -76.25, 553.0, 978.25, 491.75
    .float 587.25, -462.75, 503.5, -297.75, -835.5, -473.0, -555.75, -172.75
    .float 766.0, 199.5, -769.5, -505.5, -173.0, -86.0, 813.5, -640.75
    .float 523.5, -391.75, 309.0, 277.0, 835.75, 751.25, 842.75, 553.0
    .float 866.75, -474.25, 359.0, -841.25, 658.75, 419.0, 2.5, -128.25
    .float 720.75, 253.25, 457.5, -244.5, 745.0, -806.25, -813.75, 568.25

.align 4
twiddles:
    # For N = 8 (first N/2 twiddles)
    # Format: real0, imag0, real1, imag1, ...
    .float 1.0, 0.0      # W0
    .float 0.7071, -0.7071  # W1
    .float 0.0, -1.0        # W2
    .float -0.7071, -0.7071 # W3

## DATA DEFINE END
size: .word MatrixSize