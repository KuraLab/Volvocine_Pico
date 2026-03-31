#!/usr/bin/env python3
"""
Animate phi1, phi2 trajectory in phase space and save to MP4.
Dynamics: dphi1/dt = 2.5*pi, dphi2/dt = 2.4*pi
Phase wraps (2*pi -> 0) are shown as discontinuities in the trajectory.
"""

import numpy as np
import imageio
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from pathlib import Path
import tempfile
import os
from PIL import Image

# Parameters
omega1 = 2.5 * np.pi  # dphi1/dt
omega2 = 2.4 * np.pi  # dphi2/dt

dt = 0.01
T_end = 20.0
t = np.arange(0, T_end + dt, dt)
N = len(t)

phi1_0 = 0.0
phi2_0 = 0.0

# Export settings
export_frame_step = 4  # Sample every 4th frame
fps = 30
frame_size = (512, 512)  # width x height

output_file = Path(__file__).parent / 'phase_trajectory_phi1_phi2.mp4'

print("=== Phase Space Trajectory Animation ===")
print(f"Dynamics: phi1_dot = {omega1/np.pi:.1f}*pi, phi2_dot = {omega2/np.pi:.1f}*pi")
print(f"Time: 0 to {T_end}s (dt={dt}s, {N} steps)")
print(f"Export: every {export_frame_step}th frame, {fps} fps, size {frame_size}")

# Simulate phase evolution
print("\n[1/3] Simulating phase dynamics...")
phi1 = np.zeros(N)
phi2 = np.zeros(N)
phi1[0] = phi1_0
phi2[0] = phi2_0

for i in range(1, N):
    phi1[i] = phi1[i-1] + omega1 * dt
    phi2[i] = phi2[i-1] + omega2 * dt

# Wrap to [0, 2*pi]
phi1_wrapped = np.mod(phi1, 2 * np.pi)
phi2_wrapped = np.mod(phi2, 2 * np.pi)

# Detect phase wraps and insert NaN for discontinuities
for i in range(1, N):
    # Check if wrapped value jumped (discontinuity)
    if abs(phi1_wrapped[i] - phi1_wrapped[i-1]) > np.pi:
        phi1_wrapped[i-1] = np.nan
    if abs(phi2_wrapped[i] - phi2_wrapped[i-1]) > np.pi:
        phi2_wrapped[i-1] = np.nan

# Create list to store frames (will write all at once)
print(f"\n[2/3] Preparing video output: {output_file}")

# Animation loop
print(f"\n[3/3] Rendering and writing frames ({len(range(0, N, export_frame_step))} output frames)...")

frame_indices = np.arange(0, N, export_frame_step)
num_frames = len(frame_indices)

# Use temporary directory for frame files
import tempfile
import os
temp_dir = tempfile.mkdtemp()

for frame_num, idx in enumerate(frame_indices):
    # Create figure for each frame
    fig = Figure(figsize=(6, 6), dpi=100, frameon=False)
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(0, 2 * np.pi)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\phi_1$', fontsize=31)
    ax.set_ylabel(r'$\phi_2$', fontsize=31)
    ax.grid(True, alpha=0.3)
    
    # Plot grid lines at multiples of pi
    pi_ticks = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
    pi_labels = ['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']
    ax.set_xticks(pi_ticks)
    ax.set_xticklabels(pi_labels, fontsize=31)
    ax.set_yticks(pi_ticks)
    ax.set_yticklabels(pi_labels, fontsize=31)
    ax.tick_params(axis='both', labelsize=31)
    
    # Plot trajectory line up to current point
    indices_to_plot = np.arange(0, idx + 1)
    valid_mask = ~np.isnan(phi1_wrapped[indices_to_plot]) & ~np.isnan(phi2_wrapped[indices_to_plot])
    
    if np.any(valid_mask):
        # Plot line segments (broken at NaNs for phase wraps)
        segments_x = []
        segments_y = []
        for i in range(len(indices_to_plot)):
            if valid_mask[i]:
                segments_x.append(phi1_wrapped[indices_to_plot[i]])
                segments_y.append(phi2_wrapped[indices_to_plot[i]])
            else:
                if segments_x:
                    ax.plot(segments_x, segments_y, color='C1', linewidth=1.5, alpha=0.7)
                    segments_x = []
                    segments_y = []
        # Plot final segment
        if segments_x:
            ax.plot(segments_x, segments_y, color='C1', linewidth=1.5, alpha=0.7)
        
        # Plot current point
        ax.plot(phi1_wrapped[idx], phi2_wrapped[idx], 'ro', markersize=8)
    
    # Save frame to temporary PNG
    temp_frame = os.path.join(temp_dir, f'frame_{frame_num:05d}.png')
    fig.savefig(temp_frame, dpi=100, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    if (frame_num + 1) % 50 == 0 or frame_num + 1 == num_frames:
        progress = 100.0 * (frame_num + 1) / num_frames
        print(f"  Progress: {frame_num + 1}/{num_frames} frames ({progress:.1f}%)")

# Read all PNG frames and resize to target frame_size
print(f"\n[4/4] Loading frames and writing video...")

video_frames = []
for frame_num in range(num_frames):
    temp_frame = os.path.join(temp_dir, f'frame_{frame_num:05d}.png')
    img = Image.open(temp_frame).convert('RGB')
    img = img.resize(frame_size, Image.Resampling.LANCZOS)
    video_frames.append(np.array(img))
    os.remove(temp_frame)

os.rmdir(temp_dir)

# Write all frames to video file
print(f"\n[4/4] Writing video file using ffmpeg...")
try:
    imageio.mimsave(str(output_file), video_frames, fps=fps, codec='libx264', pixelformat='yuv420p')
    print(f"SUCCESS! Video saved: {output_file}")
except Exception as e:
    print(f"libx264 failed ({e}), trying MPEG-4 encoder...")
    try:
        imageio.mimsave(str(output_file), video_frames, fps=fps, codec='mpeg4')
    except Exception as e2:
        print(f"ERROR: Video write failed: {e2}")
        raise

print(f"SUCCESS! Video saved: {output_file}")

file_size = output_file.stat().st_size / 1024
print(f"File size: {file_size:.1f} KB")
plt.close('all')
