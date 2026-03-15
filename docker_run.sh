# enable access to xhost from the container
xhost +


# Run docker and open bash shell 
docker run -it --rm --privileged \
--env=LOCAL_USER_ID="$(id -u)" \
-v /tmp/.X11-unix:/tmp/.X11-unix:ro \
-v "/dev:/dev" \
-v $(pwd)/iiwa_qp:/catkin_ws/src/iiwa_qp:rw \
-v $(pwd)/weiss_wsg50:/catkin_ws/src/weiss_wsg50:rw \
-v $(pwd)/kuka_iiwa_support:/catkin_ws/src/kuka_iiwa_support:rw \
-v $(pwd)/launch.yml:/catkin_ws/src/launch.yml:rw \
--env="DISPLAY=$DISPLAY" \
--network host \
--name=nsp-cbf fp018/nsp-cbf-img-2:latest bash



