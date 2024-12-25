# Mohammad Abu Shams 1200549

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Define the AES S-Box.
aes_sbox = np.array([
    [0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76],
    [0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0],
    [0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15],
    [0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75],
    [0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84],
    [0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf],
    [0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8],
    [0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2],
    [0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73],
    [0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb],
    [0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79],
    [0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08],
    [0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a],
    [0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e],
    [0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf],
    [0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16]
]).flatten()


# Function to check bijection.
def check_bijection(sbox):
    return len(set(sbox)) == 256


# Walsh-Hadamard transform.
def walsh_hadamard_transform(f, n):
    N = 2 ** n
    wht = f.copy()
    for i in range(n):
        step = 2 ** i
        for j in range(0, N, step * 2):
            for k in range(step):
                a = wht[j + k]
                b = wht[j + k + step]
                wht[j + k] = a + b
                wht[j + k + step] = a - b
    return wht


# Function to calculate nonlinearity using Walsh-Hadamard transform.
def calculate_nonlinearity(sbox):
    n = 8
    min_nonlinearity = float('inf')

    for output_bit in range(n):
        f = [(x >> output_bit) & 1 for x in sbox]
        f_mapped = [1 if bit == 0 else -1 for bit in f]
        wht = walsh_hadamard_transform(f_mapped, n)
        max_wht = max(np.abs(wht))
        nonlinearity = (2 ** (n - 1)) - (0.5 * max_wht)
        min_nonlinearity = min(min_nonlinearity, nonlinearity)

    return min_nonlinearity


# Function to calculate Strict Avalanche Criterion (SAC).
def calculate_sac(sbox):
    size = 256
    input_bits = 8
    output_bits = 8
    sac_matrix = np.zeros((input_bits, output_bits))

    for i in range(size):
        for bit in range(input_bits):
            flipped_input = i ^ (1 << bit)
            output_diff = sbox[i] ^ sbox[flipped_input]
            for output_bit in range(output_bits):
                sac_matrix[bit][output_bit] += (output_diff >> output_bit) & 1

    sac_matrix /= size
    return sac_matrix


# Function to calculate Output Bits Independence Criterion (BIC).
def calculate_bic(sbox):
    size = 256
    input_bits = 8
    output_bits = 8
    bic_matrix = np.zeros((input_bits, output_bits, output_bits))

    for i in range(size):
        for bit in range(input_bits):
            flipped_input = i ^ (1 << bit)
            output_diff = sbox[i] ^ sbox[flipped_input]
            for output_bit1 in range(output_bits):
                for output_bit2 in range(output_bits):
                    if output_bit1 != output_bit2:
                        bic_matrix[bit][output_bit1][output_bit2] += (
                            ((output_diff >> output_bit1) & 1) *
                            ((output_diff >> output_bit2) & 1)
                        )

    bic_matrix /= size
    return bic_matrix


# Function to display SAC matrix as a heatmap.
def plot_sac_matrix(sac_matrix):
    plt.figure(figsize=(10, 8))
    sns.heatmap(sac_matrix, annot=True, fmt=".4f", cmap="YlGnBu")
    plt.title("Strict Avalanche Criterion (SAC) Matrix")
    plt.xlabel("Output Bits")
    plt.ylabel("Input Bits")
    plt.xticks(np.arange(8)+0.5, [f"Bit {i}" for i in range(8)], rotation=0)
    plt.yticks(np.arange(8)+0.5, [f"Bit {i}" for i in range(8)], rotation=0)
    plt.tight_layout()
    plt.show()


# Function to display BIC matrix for a specific input bit as a heatmap.
def plot_bic_matrix(bic_matrix, input_bit):
    plt.figure(figsize=(10, 8))
    sns.heatmap(bic_matrix[input_bit], annot=True, fmt=".4f", cmap="YlGnBu")
    plt.title(f"Output Bits Independence Criterion (BIC) for Input Bit {input_bit}")
    plt.xlabel("Output Bit 2")
    plt.ylabel("Output Bit 1")
    plt.xticks(np.arange(8)+0.5, [f"Bit {i}" for i in range(8)], rotation=0)
    plt.yticks(np.arange(8)+0.5, [f"Bit {i}" for i in range(8)], rotation=0)
    plt.tight_layout()
    plt.show()


# MAIN.
def main():
    print("AES S-Box Cryptographic Properties Analysis\n")

    # 1. Bijection Check.
    bijection = check_bijection(aes_sbox)
    print("1. Bijection Check:")
    print(f"   - Result: {bijection}")
    print("   - Explanation: Every input goes to only one output. This show the S-Box is bijective.\n")


    # 2. Nonlinearity Calculation.
    nonlinearity = calculate_nonlinearity(aes_sbox)
    print("2. Nonlinearity:")
    print(f"   - Minimum Nonlinearity: {nonlinearity}")
    print("   - Explanation: Big minimum nonlinearity make strong defense against linear and differential attacks.\n")


    # 3. Strict Avalanche Criterion (SAC) Calculation.
    sac = calculate_sac(aes_sbox)
    print("3. Strict Avalanche Criterion (SAC):")
    print("   - SAC Matrix:")
    # Formatting the SAC matrix for display.
    sac_table = ""
    header = "      " + "  ".join([f"Bit {i}" for i in range(8)])
    sac_table += header + "\n"
    for i, row in enumerate(sac):
        sac_table += f"Bit {i} " + "  ".join([f"{val:.4f}" for val in row]) + "\n"
    print(sac_table)
    # Plotting the SAC matrix heatmap.
    plot_sac_matrix(sac)
    print("   - Explanation: Flipping one input bit change the output bits around 50% of times.\n")


    # 4. Output Bits Independence Criterion (BIC) Calculation.
    bic = calculate_bic(aes_sbox)
    print("4. Output Bits Independence Criterion (BIC):")
    print("   - Sample BIC Results for Input Bit 0:")
    # Formatting the BIC matrix for input bit 0.
    bic_table = ""
    header = "      " + "  ".join([f"Bit {i}" for i in range(8)])
    bic_table += header + "\n"
    for i, row in enumerate(bic[0]):
        bic_table += f"Bit {i} " + "  ".join([f"{val:.4f}" for val in row]) + "\n"
    print(bic_table)
    # Plotting the BIC matrix heatmap for input bit 0
    plot_bic_matrix(bic, 0)
    print("   - Explanation: Output bits change separate when one input bit is flipped. The values need to be near 0.25.\n")


if __name__ == "__main__":
    main()
