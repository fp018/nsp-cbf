#!/usr/bin/env python3
import rospy
from geometry_msgs.msg import Pose
from gazebo_msgs.srv import SpawnModel
from tf.transformations import quaternion_from_euler

rospy.init_node("plate_spawner")

model_xml = rospy.get_param("/plate_description")

x, y, z = rospy.get_param("/plate_position")
R, P, Y = rospy.get_param("/plate_rpy", [0.0, 0.0, 0.0])

pose = Pose()
pose.position.x = x
pose.position.y = y
pose.position.z = z

qx, qy, qz, qw = quaternion_from_euler(R, P, Y)
pose.orientation.x = qx
pose.orientation.y = qy
pose.orientation.z = qz
pose.orientation.w = qw

rospy.wait_for_service("/gazebo/spawn_urdf_model")
spawn_model = rospy.ServiceProxy("/gazebo/spawn_urdf_model", SpawnModel)

resp = spawn_model(
    model_name="plate",
    model_xml=model_xml,
    robot_namespace="",
    initial_pose=pose,
    reference_frame="world"
)

rospy.loginfo(resp.status_message)
