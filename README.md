3D Ultrasound Multipass Volume Reconstruction
October 2025
April Anlage



Submitted to IEEE TMBE July 2025, revised October 2025.



5 sample phantom scans can be found here: https://www.dropbox.com/scl/fo/trtmsgp0nbmosokvof2ek/AIhzEXRbNPlq5-O1JmFQ3mA?rlkey=b0ltcvu5pfp4y2et66ewim7cw\&st=bvsyqfw5\&dl=0

Place these files in the Phantom\_data folder.



Each example includes 3 files:

(a) US\_phantom\_scan# (the DICOM ultrasound file)

(b) camPoses\_phantom\_scan# (corresponding camera poses, raw from camera)

(c) bc\_phantom\_scan# (coordinates of segmentation outlines, downsampled to US frames)





The main algorithm runs in iterations\_phantom\_r.m These file names/paths need to be edited in this file.

The rest are helper functions.

