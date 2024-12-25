# AES S-Box Cryptographic Properties Analysis


## Introduction

The Advanced Encryption Standard (AES) utilizes a substitution box (S-Box) as a fundamental component for ensuring security through confusion. This project analyzes the cryptographic properties of the AES S-Box, focusing on:

- **Bijection**: Ensuring a one-to-one mapping between input and output.
- **Nonlinearity**: Assessing resistance against linear and differential attacks.
- **Strict Avalanche Criterion (SAC)**: Measuring how input bit flipping affects output bits.
- **Output Bits Independence Criterion (BIC)**: Evaluating the independence between output bits when an input bit is flipped.

## Features

- **Bijection Check**: Validates that the AES S-Box is a bijection.
- **Nonlinearity Calculation**: Computes the minimum nonlinearity using the Walsh-Hadamard transform.
- **Strict Avalanche Criterion (SAC)**: Generates and visualizes an SAC matrix as a heatmap.
- **Output Bits Independence Criterion (BIC)**: Creates a BIC matrix for a specific input bit and visualizes it as a heatmap.

