// Author of FLOAM: Wang Han 
// Email wh200720041@gmail.com
// Homepage https://wanghan.pro
#include "laserProcessingClass.h"

void LaserProcessingClass::init(lidar::Lidar lidar_param_in){
    
    lidar_param = lidar_param_in;

}

// Lego Loam的近似实现
void LaserProcessingClass::featureExtraction(const pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_ground, const pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_not_ground, pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_out_edge, pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_out_surf){
    pcl::PointCloud<pcl::PointXYZI>::Ptr pc_temp;          
    pc_temp = pc_not_ground;
    *pc_temp += *pc_ground;
    std::sort(pc_temp->points.begin(), pc_temp->points.end(), [](const pcl::PointXYZI & a, const pcl::PointXYZI & b)
    { 
        return a.intensity < b.intensity; 
    });

    std::vector<int> indices;
    pcl::removeNaNFromPointCloud(*pc_temp, indices);
    int N_SCANS = lidar_param.num_lines;
    // 储存每个SCAN的点云 index 0是第0条SCAN的点云
    std::vector<pcl::PointCloud<pcl::PointXYZI>::Ptr> laserCloudScans;
    for(int i=0;i<N_SCANS;i++){
        // 为每圈SCAN分配空间
        laserCloudScans.push_back(pcl::PointCloud<pcl::PointXYZI>::Ptr(new pcl::PointCloud<pcl::PointXYZI>()));
    }

    // 将属于不同SCAN的点云分离开
    for (int i = 0; i < (int) pc_temp->points.size(); i++)
    {
        int scanID=0;
        // 在xy平面上的投影长度
        double distance = sqrt(pc_temp->points[i].x * pc_temp->points[i].x + pc_temp->points[i].y * pc_temp->points[i].y);
        // 滤除距离传感器太近或太远的点
        if(distance<lidar_param.min_distance || distance>lidar_param.max_distance)
            continue;
        // 和xy平面的夹角
        double angle = atan(pc_temp->points[i].z / distance) * 180 / M_PI;
        
        if (N_SCANS == 16)
        {
            // 计算该点属于第几条SCAN
            scanID = int((angle + 15) / 2 + 0.5);
            // 剔除离群点
            if (scanID > (N_SCANS - 1) || scanID < 0)
            {
                continue;
            }
        }
        else if (N_SCANS == 32)
        {
            scanID = int((angle + 92.0/3.0) * 3.0 / 4.0);
            if (scanID > (N_SCANS - 1) || scanID < 0)
            {
                continue;
            }
        }
        else if (N_SCANS == 64)
        {   
            if (angle >= -8.83)
                scanID = int((2 - angle) * 3.0 + 0.5);
            else
                scanID = N_SCANS / 2 + int((-8.83 - angle) * 2.0 + 0.5);

            if (angle > 2 || angle < -24.33 || scanID > 63 || scanID < 0)
            {
                continue;
            }
        }
        else
        {
            printf("wrong scan number\n");
        }
        laserCloudScans[scanID]->push_back(pc_temp->points[i]); 

    }

    // 对不同SCAN的点云计算它们的曲率，原理和LOAM相同
    for(int i = 0; i < N_SCANS; i++){
        if(laserCloudScans[i]->points.size()<131){
            continue;
        }
        // 一个SCAN里面所有点的曲率（除去了前后五个点的曲率）
        std::vector<Double2d> cloudCurvature; 
        int total_points = laserCloudScans[i]->points.size()-10;
        // 忽略前后五个点，对其它点计算曲率（因为前后五个点周围的点不够计算曲率）
        for(int j = 5; j < (int)laserCloudScans[i]->points.size() - 5; j++){
            // 周围有连续点
            if(fabs(laserCloudScans[i]->points[j + 5].intensity - laserCloudScans[i]->points[j - 5].intensity - 10.0) < 1e-2){
                double diffX = laserCloudScans[i]->points[j - 5].x + laserCloudScans[i]->points[j - 4].x + laserCloudScans[i]->points[j - 3].x + laserCloudScans[i]->points[j - 2].x + laserCloudScans[i]->points[j - 1].x - 10 * laserCloudScans[i]->points[j].x + laserCloudScans[i]->points[j + 1].x + laserCloudScans[i]->points[j + 2].x + laserCloudScans[i]->points[j + 3].x + laserCloudScans[i]->points[j + 4].x + laserCloudScans[i]->points[j + 5].x;
                double diffY = laserCloudScans[i]->points[j - 5].y + laserCloudScans[i]->points[j - 4].y + laserCloudScans[i]->points[j - 3].y + laserCloudScans[i]->points[j - 2].y + laserCloudScans[i]->points[j - 1].y - 10 * laserCloudScans[i]->points[j].y + laserCloudScans[i]->points[j + 1].y + laserCloudScans[i]->points[j + 2].y + laserCloudScans[i]->points[j + 3].y + laserCloudScans[i]->points[j + 4].y + laserCloudScans[i]->points[j + 5].y;
                double diffZ = laserCloudScans[i]->points[j - 5].z + laserCloudScans[i]->points[j - 4].z + laserCloudScans[i]->points[j - 3].z + laserCloudScans[i]->points[j - 2].z + laserCloudScans[i]->points[j - 1].z - 10 * laserCloudScans[i]->points[j].z + laserCloudScans[i]->points[j + 1].z + laserCloudScans[i]->points[j + 2].z + laserCloudScans[i]->points[j + 3].z + laserCloudScans[i]->points[j + 4].z + laserCloudScans[i]->points[j + 5].z;
                // 储存id和曲率
                Double2d distance(j,diffX * diffX + diffY * diffY + diffZ * diffZ);
                cloudCurvature.push_back(distance);
            }
            else{
                Double2d distance(j,-1);
                cloudCurvature.push_back(distance);
            }

        }
        // 将一个SCAN分为6个sector
        for(int j=0;j<6;j++){
            // 将一圈的点云均分为6份
            int sector_length = (int)(total_points/6);
            // 取出一个sector的点云
            int sector_start = sector_length *j;
            int sector_end = sector_length *(j+1)-1;
            // 最后一个sector拿到所有剩下的点，以防损失某些点
            if (j==5){
                sector_end = total_points - 1; 
            }
            // 一个sector的所有点
            std::vector<Double2d> subCloudCurvature(cloudCurvature.begin()+sector_start,cloudCurvature.begin()+sector_end); 
            
            featureExtractionFromSector(laserCloudScans[i],subCloudCurvature, pc_out_edge, pc_out_surf);
            
        }

    }

    pc_temp->clear();
    for(auto point : *pc_out_surf){
        if(fabs(point.intensity - int(point.intensity)) > 0.1)
            pc_temp->push_back(point);
    }

    // 将地面中的点全部作为平面点，同时加上非地面点中提取出的平面点
    // *pc_out_surf += *pc_ground;
    pc_out_surf = pc_temp;

}

// 输入：点云，每个点的曲率，线特征输出，面特征输出
void LaserProcessingClass::featureExtractionFromSector(const pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_in, std::vector<Double2d>& cloudCurvature, pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_out_edge, pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_out_surf){

    // 按照点的曲率进行排序
    std::sort(cloudCurvature.begin(), cloudCurvature.end(), [](const Double2d & a, const Double2d & b)
    { 
        return a.value < b.value; 
    });

    // 当前已经选中的点的个数
    int largestPickedNum = 0;
    // 
    std::vector<int> picked_points;
    int point_info_count =0;
    // 从曲率最大的点开始
    for (int i = cloudCurvature.size()-1; i >= 0; i--)
    {
        // 该点在pc_in中的index
        int ind = cloudCurvature[i].id; 
        // 如果该点没有被选中过
        // TODO: 专门开个数组用来记录某个点是否被选中过
        if(std::find(picked_points.begin(), picked_points.end(), ind)==picked_points.end()){
            // 曲率过小则停止搜索
            if(cloudCurvature[i].value <= 0.1){
                break;
            }
            
            largestPickedNum++;
            //标记已经选中的点
            picked_points.push_back(ind);

            // 保证一个sector中最多只有20个线特征
            if (largestPickedNum <= 20){
                pc_out_edge->push_back(pc_in->points[ind]);
                point_info_count++;
            }else{
                break;
            }
            // 滤除选中的线特征点周围的点
            for(int k=1;k<=5;k++){
                double diffX = pc_in->points[ind + k].x - pc_in->points[ind + k - 1].x;
                double diffY = pc_in->points[ind + k].y - pc_in->points[ind + k - 1].y;
                double diffZ = pc_in->points[ind + k].z - pc_in->points[ind + k - 1].z;
                // 如果两个点距离比较远，则不滤除该点
                if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05){
                    break;
                }
                picked_points.push_back(ind+k);
            }
            for(int k=-1;k>=-5;k--){
                double diffX = pc_in->points[ind + k].x - pc_in->points[ind + k + 1].x;
                double diffY = pc_in->points[ind + k].y - pc_in->points[ind + k + 1].y;
                double diffZ = pc_in->points[ind + k].z - pc_in->points[ind + k + 1].z;
                if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05){
                    break;
                }
                picked_points.push_back(ind+k);
            }

        }
    }

    //find flat points
    // point_info_count =0;
    // int smallestPickedNum = 0;
    
    // for (int i = 0; i <= (int)cloudCurvature.size()-1; i++)
    // {
    //     int ind = cloudCurvature[i].id; 

    //     if( std::find(picked_points.begin(), picked_points.end(), ind)==picked_points.end()){
    //         if(cloudCurvature[i].value > 0.1){
    //             //ROS_WARN("extracted feature not qualified, please check lidar");
    //             break;
    //         }
    //         smallestPickedNum++;
    //         picked_points.push_back(ind);
            
    //         if(smallestPickedNum <= 4){
    //             //find all points
    //             pc_surf_flat->push_back(pc_in->points[ind]);
    //             pc_surf_lessFlat->push_back(pc_in->points[ind]);
    //             point_info_count++;
    //         }
    //         else{
    //             break;
    //         }

    //         for(int k=1;k<=5;k++){
    //             double diffX = pc_in->points[ind + k].x - pc_in->points[ind + k - 1].x;
    //             double diffY = pc_in->points[ind + k].y - pc_in->points[ind + k - 1].y;
    //             double diffZ = pc_in->points[ind + k].z - pc_in->points[ind + k - 1].z;
    //             if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05){
    //                 break;
    //             }
    //             picked_points.push_back(ind+k);
    //         }
    //         for(int k=-1;k>=-5;k--){
    //             double diffX = pc_in->points[ind + k].x - pc_in->points[ind + k + 1].x;
    //             double diffY = pc_in->points[ind + k].y - pc_in->points[ind + k + 1].y;
    //             double diffZ = pc_in->points[ind + k].z - pc_in->points[ind + k + 1].z;
    //             if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05){
    //                 break;
    //             }
    //             picked_points.push_back(ind+k);
    //         }

    //     }
    // }
    
    // 从曲率最小的点开始
    for (int i = 0; i <= (int)cloudCurvature.size()-1; i++)
    {
        int ind = cloudCurvature[i].id; 
        // 如果该点没有被选中过并且该点不是无效点
        // TODO: 专门开个数组用来记录某个点是否被选中过
        if( fabs(cloudCurvature[i].value + 1) < 1e-2 && std::find(picked_points.begin(), picked_points.end(), ind)==picked_points.end() )
        {
            if(cloudCurvature[i].value > 0.05){
                break;
            }
            // 则认为该点是平面点
            pc_out_surf->push_back(pc_in->points[ind]);
        }
    }
    


}
LaserProcessingClass::LaserProcessingClass(){
    
}

Double2d::Double2d(int id_in, double value_in){
    id = id_in;
    value =value_in;
};

PointsInfo::PointsInfo(int layer_in, double time_in){
    layer = layer_in;
    time = time_in;
};
