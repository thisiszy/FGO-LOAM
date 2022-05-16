// Author of FLOAM: Wang Han 
// Email wh200720041@gmail.com
// Homepage https://wanghan.pro

//c++ lib
#include <cmath>
#include <vector>
#include <mutex>
#include <queue>
#include <thread>
#include <chrono>

//ros lib
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <nav_msgs/Odometry.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>

//pcl lib
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

//local lib
#include "lidar.h"
#include "laserProcessingClass.h"


LaserProcessingClass laserProcessing;
std::mutex mutex_lock;
std::queue<sensor_msgs::PointCloud2ConstPtr> pointCloudBuf;
lidar::Lidar lidar_param;

ros::Publisher pubEdgePoints;
ros::Publisher pubSurfPoints;
ros::Publisher pubLaserCloudFiltered;
// ros::Publisher pubGroundPoints;
// ros::Publisher pubNotGroundPoints;
PatchWork<pcl::PointXYZI> PatchworkGroundSeg;

void velodyneHandler(const sensor_msgs::PointCloud2ConstPtr &laserCloudMsg)
{
    mutex_lock.lock();
    pointCloudBuf.push(laserCloudMsg);
    mutex_lock.unlock();
   
}

double total_time =0;
int total_frame=0;

void laser_processing(){
    while(1){
        if(!pointCloudBuf.empty()){
            //read data
            mutex_lock.lock();
            pcl::PointCloud<pcl::PointXYZI>::Ptr pointcloud_in(new pcl::PointCloud<pcl::PointXYZI>());
            pcl::fromROSMsg(*pointCloudBuf.front(), *pointcloud_in);
            ros::Time pointcloud_time = (pointCloudBuf.front())->header.stamp;
            pointCloudBuf.pop();
            mutex_lock.unlock();


            std::chrono::time_point<std::chrono::system_clock> start, end;
            start = std::chrono::system_clock::now();

            // TODO: 把patchwork作为一个单独的节点发布
            // 将点云编号，编号存在intensity中
            int cloud_size = pointcloud_in->size();
            for(int i = 0; i < cloud_size; i++){
                pointcloud_in->points[i].intensity = i;
            }
            pcl::PointCloud<pcl::PointXYZI>::Ptr pointcloud_ground(new pcl::PointCloud<pcl::PointXYZI>());          
            pcl::PointCloud<pcl::PointXYZI>::Ptr pointcloud_not_ground(new pcl::PointCloud<pcl::PointXYZI>());
            double time_taken;
            PatchworkGroundSeg.estimate_ground(*pointcloud_in, *pointcloud_ground, *pointcloud_not_ground, time_taken);
            // std::sort(pointcloud_ground->points.begin(), pointcloud_ground->points.end(), [](const pcl::PointXYZI & a, const pcl::PointXYZI & b)
            // { 
            //     return a.intensity < b.intensity; 
            // });
            // 将被patchwork打乱之后的点云重新编号
            // std::sort(pointcloud_not_ground->points.begin(), pointcloud_not_ground->points.end(), [](const pcl::PointXYZI & a, const pcl::PointXYZI & b)
            // { 
            //     return a.intensity < b.intensity; 
            // });
            // 标记ground点云
            for(auto &point : pointcloud_ground->points){
                point.intensity += 0.5;
            }

            pcl::PointCloud<pcl::PointXYZI>::Ptr pointcloud_edge(new pcl::PointCloud<pcl::PointXYZI>());          
            pcl::PointCloud<pcl::PointXYZI>::Ptr pointcloud_surf(new pcl::PointCloud<pcl::PointXYZI>());
            laserProcessing.featureExtraction(pointcloud_ground, pointcloud_not_ground, pointcloud_edge, pointcloud_surf);
            
            end = std::chrono::system_clock::now();
            std::chrono::duration<float> elapsed_seconds = end - start;
            total_frame++;
            float time_temp = elapsed_seconds.count() * 1000;
            total_time+=time_temp;
            if(total_frame % 10 == 0){
                ROS_INFO("average laser processing time %f ms \n \n", total_time/total_frame);
                ROS_INFO("%d ground points, %d non-ground surface points \n \n", pointcloud_ground->size(), pointcloud_surf->size()-pointcloud_ground->size());
                ROS_INFO("%d edge points, %d surface points, total %d points \n \n", pointcloud_edge->size(), pointcloud_surf->size(), pointcloud_edge->size()+pointcloud_surf->size());
            }

            sensor_msgs::PointCloud2 laserCloudFilteredMsg;
            pcl::PointCloud<pcl::PointXYZI>::Ptr pointcloud_filtered(new pcl::PointCloud<pcl::PointXYZI>());  
            *pointcloud_filtered+=*pointcloud_edge;
            *pointcloud_filtered+=*pointcloud_surf;
            pcl::toROSMsg(*pointcloud_filtered, laserCloudFilteredMsg);
            laserCloudFilteredMsg.header.stamp = pointcloud_time;
            laserCloudFilteredMsg.header.frame_id = "base_link";
            pubLaserCloudFiltered.publish(laserCloudFilteredMsg);

            sensor_msgs::PointCloud2 edgePointsMsg;
            pcl::toROSMsg(*pointcloud_edge, edgePointsMsg);
            edgePointsMsg.header.stamp = pointcloud_time;
            edgePointsMsg.header.frame_id = "base_link";
            pubEdgePoints.publish(edgePointsMsg);


            sensor_msgs::PointCloud2 surfPointsMsg;
            pcl::toROSMsg(*pointcloud_surf, surfPointsMsg);
            surfPointsMsg.header.stamp = pointcloud_time;
            surfPointsMsg.header.frame_id = "base_link";
            pubSurfPoints.publish(surfPointsMsg);


            // sensor_msgs::PointCloud2 surfGroundMsg;
            // pcl::toROSMsg(*pointcloud_ground, surfGroundMsg);
            // surfGroundMsg.header.stamp = pointcloud_time;
            // surfGroundMsg.header.frame_id = "base_link";
            // pubGroundPoints.publish(surfGroundMsg);


            // sensor_msgs::PointCloud2 surfNotGroundMsg;
            // pcl::toROSMsg(*pointcloud_not_ground, surfNotGroundMsg);
            // surfNotGroundMsg.header.stamp = pointcloud_time;
            // surfNotGroundMsg.header.frame_id = "base_link";
            // pubNotGroundPoints.publish(surfNotGroundMsg);

        }
        //sleep 2 ms every time
        std::chrono::milliseconds dura(2);
        std::this_thread::sleep_for(dura);
    }
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "main");
    ros::NodeHandle nh;

    int scan_line = 64;
    double vertical_angle = 2.0;
    double scan_period= 0.1;
    double max_dis = 60.0;
    double min_dis = 2.0;

    string dt_file_loc;
    nh.getParam("/scan_period", scan_period); 
    nh.getParam("/vertical_angle", vertical_angle); 
    nh.getParam("/max_dis", max_dis);
    nh.getParam("/min_dis", min_dis);
    nh.getParam("/scan_line", scan_line);
    nh.getParam("/dt_folder", dt_file_loc);

    lidar_param.setScanPeriod(scan_period);
    lidar_param.setVerticalAngle(vertical_angle);
    lidar_param.setLines(scan_line);
    lidar_param.setMaxDistance(max_dis);
    lidar_param.setMinDistance(min_dis);

    laserProcessing.init(lidar_param);
    PatchworkGroundSeg.init(nh);

    ros::Subscriber subLaserCloud = nh.subscribe<sensor_msgs::PointCloud2>("/velodyne_points", 100, velodyneHandler);

    pubLaserCloudFiltered = nh.advertise<sensor_msgs::PointCloud2>("/velodyne_points_filtered", 100);

    pubEdgePoints = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_edge", 100);

    pubSurfPoints = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_surf", 100); 

    // pubGroundPoints = nh.advertise<sensor_msgs::PointCloud2>("/laser_ground", 100); 

    // pubNotGroundPoints = nh.advertise<sensor_msgs::PointCloud2>("/laser_not_ground", 100); 

    std::thread laser_processing_process{laser_processing};

    ROS_INFO("\033[1;32m---->\033[0m Laser Processing Started.");

    ros::spin();

    FILE * time_file = fopen((dt_file_loc+"time_proc.txt").c_str(), "w");
    fprintf(time_file, "%f\n", total_time/total_frame);
    fclose(time_file);

    printf("\033[1;33maverage laser processing time %f ms\033[0m \n \n", total_time/total_frame);

    return 0;
}

