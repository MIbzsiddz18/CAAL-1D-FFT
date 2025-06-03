import sys
import math
import random

def generate_test_signal(size, signal_type="cosine"):
    """Generate different types of test signals for FFT"""
    
    if signal_type == "cosine":
        # Single frequency cosine wave
        frequency = 1  # One cycle over the signal length
        signal = []
        for i in range(size):
            angle = 2 * math.pi * frequency * i / size
            signal.append(math.cos(angle))
        return signal
    
    elif signal_type == "sine":
        # Single frequency sine wave
        frequency = 1
        signal = []
        for i in range(size):
            angle = 2 * math.pi * frequency * i / size
            signal.append(math.sin(angle))
        return signal
    
    elif signal_type == "multi_sine":
        # Multiple frequency sine waves
        frequencies = [1, 3, 5]  # Multiple harmonics
        amplitudes = [1.0, 0.5, 0.25]
        signal = []
        for i in range(size):
            value = 0
            for freq, amp in zip(frequencies, amplitudes):
                angle = 2 * math.pi * freq * i / size
                value += amp * math.sin(angle)
            signal.append(value)
        return signal
    
    elif signal_type == "impulse":
        # Impulse signal (delta function)
        signal = [0.0] * size
        signal[0] = 1.0  # Impulse at the beginning
        return signal
    
    elif signal_type == "step":
        # Step function
        signal = []
        for i in range(size):
            if i < size // 2:
                signal.append(0.0)
            else:
                signal.append(1.0)
        return signal
    
    elif signal_type == "square":
        # Square wave
        frequency = 1
        signal = []
        for i in range(size):
            angle = 2 * math.pi * frequency * i / size
            signal.append(1.0 if math.sin(angle) >= 0 else -1.0)
        return signal
    
    elif signal_type == "noise":
        # White noise
        signal = []
        for i in range(size):
            signal.append(random.uniform(-1.0, 1.0))
        return signal
    
    elif signal_type == "chirp":
        # Frequency sweep (chirp signal)
        f_start = 0.1
        f_end = 5.0
        signal = []
        for i in range(size):
            t = i / size
            freq = f_start + (f_end - f_start) * t
            angle = 2 * math.pi * freq * t * size
            signal.append(math.cos(angle))
        return signal
    
    else:
        raise ValueError(f"Unknown signal type: {signal_type}")

def print_signal_info(signal, signal_type):
    """Print basic information about the generated signal"""
    print(f"Signal type: {signal_type}")
    print(f"Length: {len(signal)}")
    print(f"Min value: {min(signal):.4f}")
    print(f"Max value: {max(signal):.4f}")
    print(f"Mean value: {sum(signal)/len(signal):.4f}")
    print()

def main():
    """Test the signal generation with different types"""
    size = 64  # Signal length
    
    signal_types = ["cosine", "sine", "multi_sine", "impulse", 
                   "step", "square", "noise", "chirp"]
    
    for sig_type in signal_types:
        signal = generate_test_signal(size, sig_type)
        print_signal_info(signal, sig_type)
        
        # Print first few samples for inspection
        print(f"First 8 samples: {[f'{x:.4f}' for x in signal[:8]]}")
        print("-" * 50)

if __name__ == "__main__":
    main()
