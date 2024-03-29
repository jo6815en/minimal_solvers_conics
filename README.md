# Fast and efficient minimal solvers for quadric based camera pose estimation

*Abstract* In this paper we address absolute camera pose estimation. An efficient (and standard) way to solve this problem, is to use sparse keypoint correspondences. In many cases point features are not available, or are unstable over time and  viewing conditions.  We propose a framework based on silhouettes of quadric surfaces, with special emphasis on cylinders.  We provide mathematical analysis of the problem of projected cylinders in particular, but also general quadrics. We develop a number of minimal solvers for estimating camera pose from silhouette lines of cylinders, given different calibration and cylinder properties. These solvers can be used efficiently in  bootstrapping robust estimation schemes, such as RANSAC.  Note that even though we have lines as image features, this is a different case than line based pose estimation, since we do not have 2D-line to 3D-line correspondences. We perform synthetic accuracy and robustness tests and evaluate on a number of real case scenarios. 

## Gererate examples

In MATLAB:

```bash
# Random oriented cylinders
[Roterr, transerr] = conics_synt
```

The script ``conics_synt.m`` simulate random conics and returns the rotation and translation error based on
the estimated camera matrix compared to the true one.

![image](images/synt_image.jpg)

```bash
# Parallel cylinders
[Roterr, transerr] = conics_synt_parallel
```

The script ``conics_synt_parallel.m`` simulate random parallel conics and returns the rotation and translation error based on
the estimated camera matrix compared to the true one.

![image](images/synt_parallell_image.jpg)

```bash
# Unknown focal length of the camera
roller_coaster_pose
```

The script ``roller_coaster_pose.m`` takes in an images of a roller coaster captured with a camera without known focal length and returns a camera matrix, focal length and images showing the reprojected lines of the silhouettes of the cylinders. 

![image](images/rollercoaster.jpg)

## Cite
```bash
@inproceedings{gummesonengman2022fast,
  title={Fast and efficient minimal solvers for quadric based camera pose estimation},
  author={Gummeson, Anna and Engman, Johanna and {\AA}str{\"o}m, Kalle and Oskarsson, Magnus},
  booktitle={Proceedings of the International Conference on Pattern Recognition},
  year={2022}
}
```
