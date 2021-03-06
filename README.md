# Relative-Pose-Estimation-for-LF-Cameras
Relative Pose Estimation for Light Field Cameras Based on LF-Point-LF-Point Correspondence Model

## Preparation

- VLFeat ([Download link](https://www.vlfeat.org/download.html))

Please firstly setup VLFeat in MATLAB (at least 2009B) successfully, and then try to run the codes. The instructions are presented [here](https://www.vlfeat.org/install-matlab.html).

- LFs ([Download link](https://pan.baidu.com/s/1qfivKb8pYvaIGuGR8NvGPg)) [Extraction code: LFPE]

Please download files "IMG_1__Decoded.mat" and "IMG_2__Decoded.mat", and put these two files in the folder "data\monkeyking\". The folder structure is as

```tex
data/
├──monkeyking/
          ├──depth/…
          ├──IMG_1__Decoded.mat
          ├──IMG_1__Decoded.mat
```

Note that two files have the same name "IMG_1__Decoded.mat" (also two files have the same name "IMG_2__Decoded.mat"), one is in the folder "data\monkeyking\", and the other is in the folder "data\monkeyking\depth\E4\". Do not mix up them.

## Run codes

Run Main.m

## Performance

![](https://github.com/SaipingZhang/Relative-Pose-Estimation-for-LF-Cameras/blob/main/performance/monkeyking.png)


If you have any question or find any bug, please feel free to contact:

Saiping Zhang

Email: spzhang@stu.xidian.edu.cn
