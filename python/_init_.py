import math
import struct
import sys
import matplotlib.pyplot as plt
import numpy as np

def find_and_save_complex_lines(filename):
    """Extract complex number data from log file"""
    all_results = []
    result = []
    size_line = ""
    is_saving = False

    with open(filename, 'r') as file:
        prev_line = None

        for line in file:
            line_strip = line.strip()
            columns = line_strip.split()

            if prev_line:
                prev_columns = prev_line.strip().split()

                if len(prev_columns) > 6 and len(columns) > 6:
                    if prev_columns[6] == "00000123" and columns[6] == "00000456":
                        if is_saving:
                            all_results.append(result)
                            result = []
                            is_saving = False
                        else:
                            is_saving = True

            if is_saving:
                if "c.mv     a1" in line_strip:
                    size_line = line_strip
                result.append(line_strip)

            prev_line = line

    if result:
        all_results.append(result)

    return all_results, size_line

def filter_lines_by_flw(lines):
    """Filter lines containing floating point load instructions"""
    return [line for line in lines if "flw" in line]

def extract_7th_column(lines):
    """Extract the 7th column which contains hex values"""
    return [line.strip().split()[6] for line in lines if len(line.strip().split()) > 6]

def hex_to_float(hex_array):
    """Convert hex strings to float values"""
    result = []
    
    for hex_str in hex_array:
        try:
            hex_value = int(hex_str, 16)
            float_value = struct.unpack('!f', struct.pack('!I', hex_value))[0]
            result.append(float_value)
        except ValueError:
            pass
    
    return result

def parse_complex_array(float_array, size):
    """Convert flat float array to complex number array"""
    complex_array = []
    for i in range(0, len(float_array), 2):
        if i + 1 < len(float_array):
            real = float_array[i]
            imag = float_array[i + 1]
            complex_array.append(complex(real, imag))
        else:
            complex_array.append(complex(float_array[i], 0))
    
    return complex_array[:size]

def compute_reference_fft(signal):
    """Compute reference FFT using NumPy for comparison"""
    return np.fft.fft(signal)

def print_complex_arrays(arrays, size, labels=None):
    """Print complex arrays in a readable format"""
    if labels is None:
        labels = [f"Array {i+1}" for i in range(len(arrays))]
    
    for i, array in enumerate(arrays):
        print(f"\n{labels[i]}:")
        print("-" * 50)
        for j, val in enumerate(array):
            magnitude = abs(val)
            phase = math.atan2(val.imag, val.real) * 180 / math.pi
            print(f"[{j:2d}] {val.real:10.6f} + {val.imag:10.6f}i  "
                  f"(mag: {magnitude:8.6f}, phase: {phase:7.2f}°)")

def plot_fft_results(input_signal, output_signal, reference_fft, size):
    """Create visualization plots for FFT results"""
    
    # Create frequency bins
    freqs = np.fft.fftfreq(size, 1.0)
    
    # Create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # Plot 1: Input signal (time domain)
    time_indices = range(size)
    real_parts = [x.real for x in input_signal]
    imag_parts = [x.imag for x in input_signal]
    
    ax1.stem(time_indices, real_parts, linefmt='b-', markerfmt='bo', label='Real')
    ax1.stem(time_indices, imag_parts, linefmt='r-', markerfmt='ro', label='Imaginary')
    ax1.set_title('Input Signal (Time Domain)')
    ax1.set_xlabel('Sample Index')
    ax1.set_ylabel('Amplitude')
    ax1.legend()
    ax1.grid(True)
    
    # Plot 2: FFT Output Magnitude
    output_magnitude = [abs(x) for x in output_signal]
    reference_magnitude = [abs(x) for x in reference_fft]
    
    ax2.stem(freqs[:size//2], output_magnitude[:size//2], 
             linefmt='b-', markerfmt='bo', label='Assembly FFT')
    ax2.stem(freqs[:size//2], reference_magnitude[:size//2], 
             linefmt='r--', markerfmt='ro', label='NumPy FFT', alpha=0.7)
    ax2.set_title('FFT Magnitude Spectrum')
    ax2.set_xlabel('Frequency Bin')
    ax2.set_ylabel('Magnitude')
    ax2.legend()
    ax2.grid(True)
    
    # Plot 3: FFT Output Phase
    output_phase = [math.atan2(x.imag, x.real) * 180 / math.pi for x in output_signal]
    reference_phase = [math.atan2(x.imag, x.real) * 180 / math.pi for x in reference_fft]
    
    ax3.stem(freqs[:size//2], output_phase[:size//2], 
             linefmt='b-', markerfmt='bo', label='Assembly FFT')
    ax3.stem(freqs[:size//2], reference_phase[:size//2], 
             linefmt='r--', markerfmt='ro', label='NumPy FFT', alpha=0.7)
    ax3.set_title('FFT Phase Spectrum')
    ax3.set_xlabel('Frequency Bin')
    ax3.set_ylabel('Phase (degrees)')
    ax3.legend()
    ax3.grid(True)
    
    # Plot 4: Error Analysis
    magnitude_error = [abs(a - b) for a, b in zip(output_magnitude, reference_magnitude)]
    ax4.stem(range(size), magnitude_error, linefmt='g-', markerfmt='go')
    ax4.set_title('Magnitude Error (Assembly vs NumPy)')
    ax4.set_xlabel('Frequency Bin')
    ax4.set_ylabel('Absolute Error')
    ax4.grid(True)
    
    plt.tight_layout()
    plt.show()

def analyze_fft_accuracy(assembly_fft, reference_fft):
    """Analyze the accuracy of the assembly FFT implementation"""
    
    print("\n" + "="*60)
    print("FFT ACCURACY ANALYSIS")
    print("="*60)
    
    # Calculate errors
    magnitude_errors = []
    phase_errors = []
    
    for i, (asm_val, ref_val) in enumerate(zip(assembly_fft, reference_fft)):
        mag_error = abs(abs(asm_val) - abs(ref_val))
        magnitude_errors.append(mag_error)
        
        asm_phase = math.atan2(asm_val.imag, asm_val.real)
        ref_phase = math.atan2(ref_val.imag, ref_val.real)
        phase_error = abs(asm_phase - ref_phase)
        if phase_error > math.pi:
            phase_error = 2 * math.pi - phase_error
        phase_errors.append(phase_error * 180 / math.pi)
    
    # Statistics
    max_mag_error = max(magnitude_errors)
    avg_mag_error = sum(magnitude_errors) / len(magnitude_errors)
    max_phase_error = max(phase_errors)
    avg_phase_error = sum(phase_errors) / len(phase_errors)
    
    print(f"Maximum magnitude error: {max_mag_error:.6f}")
    print(f"Average magnitude error: {avg_mag_error:.6f}")
    print(f"Maximum phase error: {max_phase_error:.2f}°")
    print(f"Average phase error: {avg_phase_error:.2f}°")
    
    # Check if errors are within acceptable range
    mag_threshold = 1e-3
    phase_threshold = 5.0  # degrees
    
    if max_mag_error < mag_threshold and max_phase_error < phase_threshold:
        print("\n✓ FFT implementation appears to be CORRECT!")
    else:
        print("\n✗ FFT implementation may have issues.")
        print(f"  Magnitude threshold: {mag_threshold}")
        print(f"  Phase threshold: {phase_threshold}°")

def main():
    if len(sys.argv) != 3:
        print("Usage: python fft_visualizer.py <size> <type>")
        print("Example: python fft_visualizer.py 8 NV")
        sys.exit(1)
    
    size = int(sys.argv[1])
    type_arg = sys.argv[2]
    
    # Determine log file based on type
    if type_arg == "NV":
        filename = "veer/tempFiles/logNV.txt"
    elif type_arg == "V":
        filename = "veer/tempFiles/logV.txt"
    else:
        print("Type must be 'NV' (non-vectorized) or 'V' (vectorized)")
        sys.exit(1)
    
    try:
        # Extract data from log file
        arrays, size_line = find_and_save_complex_lines(filename)
        
        if not arrays or not size_line:
            print("No FFT data found in log file.")
            sys.exit(1)
        
        # Extract size from log
        extracted_size = int(extract_7th_column([size_line])[0], 16)
        print(f"FFT size detected: {extracted_size}")
        
        # Process the arrays
        complex_arrays = []
        labels = ["Input Signal", "FFT Output"]
        
        for i, array in enumerate(arrays):
            float_data = hex_to_float(extract_7th_column(filter_lines_by_flw(array)))
            complex_data = parse_complex_array(float_data, extracted_size)
            complex_arrays.append(complex_data)
        
        if len(complex_arrays) >= 2:
            input_signal = complex_arrays[0]
            output_signal = complex_arrays[1]
            
            # Compute reference FFT
            reference_fft = compute_reference_fft(input_signal)
            
            # Print results
            print_complex_arrays([input_signal, output_signal, reference_fft], 
                                extracted_size, 
                                ["Input Signal", "Assembly FFT Output", "NumPy Reference FFT"])
            
            # Analyze accuracy
            analyze_fft_accuracy(output_signal, reference_fft)
            
            # Create visualization
            try:
                plot_fft_results(input_signal, output_signal, reference_fft, extracted_size)
            except Exception as e:
                print(f"Plotting failed: {e}")
                print("Install matplotlib and numpy for visualization: pip install matplotlib numpy")
        
        else:
            print("Insufficient data arrays found. Expected at least 2 (input and output).")
    
    except FileNotFoundError:
        print(f"Log file not found: {filename}")
        print("Make sure to run the FFT program first.")
    except Exception as e:
        print(f"Error processing data: {e}")

if __name__ == "__main__":
    main()
