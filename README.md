# nsp-cbf
This repository contains the code for the **Gazebo-ROS-Noetic** simulation of an **LBR iiwa** manipulator and the execution of a time-varying prioritized stack of tasks. The QP problem is solved with the library *osqp* [1], while the Aruco detection uses *visp* [2].

# To Simulate

A docker image is available to ease the code running. To simulate, run

```bash
# Run the container. The command will also pull the docker image from DockerHub
./docker_run.sh
```
```bash
# Compile the workspace
cd /catkin_ws && catkin build
```
```bash
# Open multiple terminals and start the simulation
tmuxp load launch.yml
```
To exit from the terminal multiplexer
```bash
# In one terminal
tmux kill-session
```


# References

[1] B. Stellato, G. Banjac, P. Goulart, A. Bemporad, and S. Boyd, “OSQP: an operator splitting solver for quadratic programs,” Mathematical Programming Computation, vol. 12, no. 4, pp. 637–672, 2020, doi: 10.1007/s12532-020-00179-2.

[2] E. Marchand, F. Spindler, F. Chaumette. ViSP for visual servoing: a generic software platform with a wide class of robot control skills. IEEE Robotics and Automation Magazine, Special Issue on “Software Packages for Vision-Based Control of Motion”, P. Oh, D. Burschka (Eds.), 12(4):40-52, December 2005.
