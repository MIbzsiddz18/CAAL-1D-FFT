/**
 * 1D Fast Fourier Transform (FFT) Implementation
 * Based on the Radix-2 Decimation-In-Time (DIT) algorithm
 * 
 * This implementation follows the matrix-based approach described in the provided documentation.
 * It computes the FFT of a complex-valued signal of length N, where N is a power of 2.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846

/**
 * Structure to represent a complex number
 * Can be used instead of complex.h if not available
 */
typedef struct {
    double real;
    double imag;
} Complex;

/**
 * Perform bit reversal permutation on the input array
 * 
 * @param x Input/output array
 * @param N Size of the array (must be a power of 2)
 */
void bit_reverse_permutation(Complex *x, int N) {
    int i, j, k;
    Complex temp;
    
    // Bit reversal
    j = 0;
    for (i = 0; i < N - 1; i++) {
        if (i < j) {
            // Swap x[i] and x[j]
            temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
        k = N / 2;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

/**
 * Compute the complex multiplication: a * b
 * 
 * @param a First complex number
 * @param b Second complex number
 * @return Result of a * b
 */
Complex complex_multiply(Complex a, Complex b) {
    Complex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

/**
 * Compute the FFT of a complex-valued signal
 * 
 * @param x Input/output array (in-place computation)
 * @param N Size of the array (must be a power of 2)
 */
void fft(Complex *x, int N) {
    int stage, butterfly, group, pair;
    int butterflySize, groupSize, distance;
    double angle;
    Complex twiddle, temp;
    
    // Perform bit reversal permutation
    bit_reverse_permutation(x, N);
    
    // FFT computation using butterfly operations
    butterflySize = 2;  // Start with 2-point DFT
    
    // Loop for log2(N) stages
    for (stage = 0; stage < log2(N); stage++) {
        groupSize = butterflySize / 2;
        distance = butterflySize;
        
        // Loop through groups
        for (group = 0; group < groupSize; group++) {
            // Compute twiddle factor W_N^group
            angle = -2.0 * PI * group / butterflySize;
            twiddle.real = cos(angle);
            twiddle.imag = sin(angle);
            
            // Loop through butterflies
            for (butterfly = 0; butterfly < N; butterfly += distance) {
                pair = butterfly + groupSize;
                
                // Butterfly operation
                temp = complex_multiply(x[pair], twiddle);
                x[pair].real = x[butterfly].real - temp.real;
                x[pair].imag = x[butterfly].imag - temp.imag;
                x[butterfly].real += temp.real;
                x[butterfly].imag += temp.imag;
            }
        }
        
        // Double butterfly size for next stage
        butterflySize *= 2;
    }
}

/**
 * Print complex array
 * 
 * @param x Array to print
 * @param N Size of the array
 * @param label Label for the output
 */
void print_complex_array(Complex *x, int N, const char *label) {
    int i;
    printf("%s:\n", label);
    for (i = 0; i < N; i++) {
        printf("[%d] %.4f + %.4fi\n", i, x[i].real, x[i].imag);
    }
    printf("\n");
}

/**
 * Generate a cosine wave as a test signal
 * 
 * @param x Output array
 * @param N Size of the array
 * @param frequency Frequency of the cosine wave
 * @param sampleRate Sample rate
 */
void generate_cosine_wave(Complex *x, int N, double frequency, double sampleRate) {
    int i;
    double t;
    
    for (i = 0; i < N; i++) {
        t = (double)i / sampleRate;
        x[i].real = cos(2.0 * PI * frequency * t);
        x[i].imag = 0.0;  // Real signal
    }
}

/**
 * Generate a complex test signal with both real and imaginary parts
 * 
 * @param x Output array
 * @param N Size of the array
 * @param frequency1 Frequency for real part
 * @param frequency2 Frequency for imaginary part
 * @param sampleRate Sample rate
 */
void generate_complex_signal(Complex *x, int N, double frequency1, double frequency2, double sampleRate) {
    int i;
    double t;
    
    for (i = 0; i < N; i++) {
        t = (double)i / sampleRate;
        x[i].real = cos(2.0 * PI * frequency1 * t);
        x[i].imag = sin(2.0 * PI * frequency2 * t);
    }
}

/**
 * Check if a number is a power of 2
 * 
 * @param n Number to check
 * @return 1 if n is a power of 2, 0 otherwise
 */
int is_power_of_two(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

/**
 * Main function to demonstrate FFT
 */
int main() {
    int N = 1024;  // Signal length (must be a power of 2)
    double sampleRate = 1024.0;  // Sample rate in Hz
    double frequency = 100.0;    // Frequency of test signal (100 Hz as mentioned in the document)
    Complex *signal;
    
    // Check if N is a power of 2
    if (!is_power_of_two(N)) {
        printf("Error: Signal length must be a power of 2.\n");
        return 1;
    }
    
    // Allocate memory for the signal
    signal = (Complex *)malloc(N * sizeof(Complex));
    if (signal == NULL) {
        printf("Error: Memory allocation failed.\n");
        return 1;
    }
    
    // Test 1: Generate and process a real cosine wave
    printf("TEST 1: REAL COSINE WAVE\n");
    printf("Generating a %d-point cosine wave with frequency %.1f Hz\n", N, frequency);
    generate_cosine_wave(signal, N, frequency, sampleRate);
    
    // Print first 10 samples of the input signal
    printf("First 10 samples of the input signal:\n");
    for (int i = 0; i < 10; i++) {
        printf("[%d] %.4f + %.4fi\n", i, signal[i].real, signal[i].imag);
    }
    
    // Compute FFT
    printf("Computing FFT...\n");
    fft(signal, N);
    
    // Print first 10 samples of the FFT result
    printf("First 10 samples of the FFT result:\n");
    for (int i = 0; i < 10; i++) {
        printf("[%d] %.4f + %.4fi (magnitude: %.4f)\n", 
               i, signal[i].real, signal[i].imag, 
               sqrt(signal[i].real * signal[i].real + signal[i].imag * signal[i].imag));
    }
    
    // Test 2: Generate and process a complex signal
    printf("\nTEST 2: COMPLEX SIGNAL\n");
    printf("Generating a %d-point complex signal\n", N);
    generate_complex_signal(signal, N, frequency, frequency/2, sampleRate);
    
    // Print first 10 samples of the input signal
    printf("First 10 samples of the input signal:\n");
    for (int i = 0; i < 10; i++) {
        printf("[%d] %.4f + %.4fi\n", i, signal[i].real, signal[i].imag);
    }
    
    // Compute FFT
    printf("Computing FFT...\n");
    fft(signal, N);
    
    // Print first 10 samples of the FFT result
    printf("First 10 samples of the FFT result:\n");
    for (int i = 0; i < 10; i++) {
        printf("[%d] %.4f + %.4fi (magnitude: %.4f)\n", 
               i, signal[i].real, signal[i].imag, 
               sqrt(signal[i].real * signal[i].real + signal[i].imag * signal[i].imag));
    }
    
    // Free memory
    free(signal);
    
    return 0;
}
