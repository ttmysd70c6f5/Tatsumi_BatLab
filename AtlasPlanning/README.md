# AtlasPlanning
This repository is for planning and visualizing a probe trajectory. `AtlasTrajectoryV5.ipynb` is the main script, and the user needs to run only this script.

## Required python packages
- opencv: run `pip install opencv-python` to install
- Pillow: run `pip install --upgrade Pillow` to install
- Matplotlib: run `conda install matplotlib` to install

## How to use
First of all, you need to import the `TrajectoryPlanning_v2` function from `TrajectoryPlanning.py`.
```python
from TrajectoryPlanning import TrajectoryPlanning_v2 as trj
```
To run the function, you need to specify two target regions from the atlas. You need to give the plate number (AP) and the coordinates (ML and DV) of the target regions in the atlas. The images of the atlas are contained in the **atlas** directory.

In the example below, the hippocampal CA1 (Plate = 52, ML = 3.0 mm, DV = 3.0 mm) and the MEC (Plate = 76, ML = 6.4 mm, DV = 10.0 mm) are targeted.

```python
# 1st target region
plate1    = 52 # plate number
ML1_mm = 3.0 # ML coordinate 
DV1_mm = 3.0 # DV coordinate
# 2nd target region
plate2    = 76 # plate number 
ML2_mm = 6.4 # ML coordinate
DV2_mm = 10.0 # DV coordinate
```
![0052](https://github.com/ttmysd70c6f5/Tatsumi_BatLab/assets/61156941/536a689b-4862-4d13-aea3-36dc9278f283)
![0076](https://github.com/ttmysd70c6f5/Tatsumi_BatLab/assets/61156941/7a4ea1f5-85b5-4770-9bbb-3c8b87febb41)

You also need to name the trajectory. The specified name is used for the name of the output directory.
```python
trajectory_name = 'Example_trajectory' # name of output directory
```

Then, run the function to compute and visualize the trajectory.
```python
trj(trajectory_name, plate1, plate2, ML1_mm, ML2_mm, DV1_mm, DV2_mm)
```

The information about the target regions and the probe angle (`Azimuth` and `Pitch`) is output. 
```python
Target region 1: Plate 81 , sinus 1.28 (mm)
Target region 2: Plate 84 , sinus 1.03 (mm)
Azimuth: 83.05 degree
Pitch: 17.01 degree
```

A new directory with the same name as `trajectory_name` is made within the working directory. This output directory contains new atlas images with the plot of the estimated probe position as a blue circle. No circle is plotted on the plate that the estimated position is out of range (`ML = [0,8]`, `DV = [0,11.6]`). 
The probe angle and the trajectory length from the 1st plate are plotted at the top of each plate.
![054_AP5 38](https://github.com/ttmysd70c6f5/Tatsumi_BatLab/assets/61156941/660409c2-6108-43cb-91e8-6ae44c5d3f7b)

If both of the two target regions are on the same plate, only that plate is processed and output in the output directory. Please note that the **length** shown in the output image is the distance between two target regions in this mode (c.f. the trajectory length from the first plate is shown if the target regions are on different plates).

```python
from TrajectoryPlanning import TrajectoryPlanning_v2 as trj
# 1st target region
plate1    = 81 # plate number
ML1_mm = 5.15 # ML coordinate 
DV1_mm = 2.85 # DV coordinate
# 2nd target region
plate2    = 81 # plate number 
ML2_mm = 7.2 # ML coordinate
DV2_mm = 9.6 # DV coordinate

trajectory_name = 'Example_trajectory2' # name of output directory

trj(trajectory_name, plate1, plate2, ML1_mm, ML2_mm, DV1_mm, DV2_mm) # run the function
```
![081_AP1 28](https://github.com/ttmysd70c6f5/Tatsumi_BatLab/assets/61156941/03c98002-6b7f-47f5-bca6-d22bfc70cdb2)
