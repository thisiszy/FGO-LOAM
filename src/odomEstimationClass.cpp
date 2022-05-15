// Author of FLOAM: Wang Han 
// Email wh200720041@gmail.com
// Homepage https://wanghan.pro

#include "odomEstimationClass.h"

void OdomEstimationClass::init(lidar::Lidar lidar_param, double map_resolution){
    //init local map
    laserCloudCornerMap = pcl::PointCloud<pcl::PointXYZI>::Ptr(new pcl::PointCloud<pcl::PointXYZI>());
    laserCloudSurfMap = pcl::PointCloud<pcl::PointXYZI>::Ptr(new pcl::PointCloud<pcl::PointXYZI>());

    //downsampling size
    downSizeFilterEdge.setLeafSize(map_resolution, map_resolution, map_resolution);
    downSizeFilterSurf.setLeafSize(map_resolution * 2, map_resolution * 2, map_resolution * 2);

    //kd-tree
    #ifdef NANOFLANN
    #else
    kdtreeEdgeMap = pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr(new pcl::KdTreeFLANN<pcl::PointXYZI>());
    kdtreeSurfMap = pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr(new pcl::KdTreeFLANN<pcl::PointXYZI>());
    #endif

    odom = Eigen::Isometry3d::Identity();
    last_odom = Eigen::Isometry3d::Identity();
    optimization_count=2;
}

void OdomEstimationClass::initMapWithPoints(const pcl::PointCloud<pcl::PointXYZI>::Ptr& edge_in, const pcl::PointCloud<pcl::PointXYZI>::Ptr& surf_in){
    // 添加特征点到Map
    *laserCloudCornerMap += *edge_in;
    *laserCloudSurfMap += *surf_in;
    optimization_count=12;
}


void OdomEstimationClass::updatePointsToMap(const pcl::PointCloud<pcl::PointXYZI>::Ptr& edge_in, const pcl::PointCloud<pcl::PointXYZI>::Ptr& surf_in){

    // 为什么需要把优化次数慢慢降下来？？
    if(optimization_count>2)
        optimization_count--;

    // 利用前一次位姿变换来预测当前位姿变换 odom * last_odom.inverse 相当于odom相对于last_odom的增量
    Eigen::Isometry3d odom_prediction = odom * (last_odom.inverse() * odom);
    last_odom = odom;
    odom = odom_prediction;

    // 储存当前点云的位姿
    q_w_curr = Eigen::Quaterniond(odom.rotation());
    t_w_curr = odom.translation();

    pcl::PointCloud<pcl::PointXYZI>::Ptr downsampledEdgeCloud(new pcl::PointCloud<pcl::PointXYZI>());
    pcl::PointCloud<pcl::PointXYZI>::Ptr downsampledSurfCloud(new pcl::PointCloud<pcl::PointXYZI>());
    // 下采样两个特征点云
    downSamplingToMap(edge_in,downsampledEdgeCloud,surf_in,downsampledSurfCloud);
    //ROS_WARN("point nyum%d,%d",(int)downsampledEdgeCloud->points.size(), (int)downsampledSurfCloud->points.size());
    if(laserCloudCornerMap->points.size()>10 && laserCloudSurfMap->points.size()>50){
        #ifdef NANOFLANN
        kdtreeEdgeMap.reset();
        kdtreeEdgeMap = std::make_unique<PC2KD>(laserCloudCornerMap);
        kdtreeSurfMap.reset();
        kdtreeSurfMap = std::make_unique<PC2KD>(laserCloudSurfMap);
        #else
        kdtreeEdgeMap->setInputCloud(laserCloudCornerMap);
        kdtreeSurfMap->setInputCloud(laserCloudSurfMap);
        #endif

        // 优化optimization_count轮（每次会对parameter在前一次基础上进行优化）
        for (int iterCount = 0; iterCount < optimization_count; iterCount++){
            {
                ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
                ceres::Problem::Options problem_options;
                // 优化问题
                ceres::Problem problem(problem_options);

                // 设置需要优化的参数
                problem.AddParameterBlock(paramEuler, 6);
                
                // 添加面特征项
                addSurfCostFactor(downsampledSurfCloud,laserCloudSurfMap,problem,loss_function);

                // 优化参数设置
                ceres::Solver::Options options;
                options.linear_solver_type = ceres::DENSE_QR;
                options.max_num_iterations = 4;
                options.minimizer_progress_to_stdout = false;
                options.check_gradients = false;
                options.gradient_check_relative_precision = 1e-4;
                ceres::Solver::Summary summary;

                // solve 会修改paramEuler
                ceres::Solve(options, &problem, &summary);
                updatePose();
            }
        }
        for (int iterCount = 0; iterCount < optimization_count; iterCount++){
            {
                ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
                ceres::Problem::Options problem_options;
                // 优化问题
                ceres::Problem problem(problem_options);

                // 设置需要优化的参数
                problem.AddParameterBlock(paramEuler, 6);
                
                // 添加线特征项
                addEdgeCostFactor(downsampledEdgeCloud,laserCloudCornerMap,problem,loss_function);

                // 优化参数设置
                ceres::Solver::Options options;
                options.linear_solver_type = ceres::DENSE_QR;
                options.max_num_iterations = 4;
                options.minimizer_progress_to_stdout = false;
                options.check_gradients = false;
                options.gradient_check_relative_precision = 1e-4;
                ceres::Solver::Summary summary;

                // solve 会修改paramEuler
                ceres::Solve(options, &problem, &summary);
                updatePose();
            }

        }
    }else{
        printf("not enough points in map to associate, map error");
    }
    // 更新odom坐标，为下一次prediction作准备
    odom = Eigen::Isometry3d::Identity();
    odom.linear() = q_w_curr.toRotationMatrix();
    odom.translation() = t_w_curr;
    // 添加特征点到map中
    addPointsToMap(downsampledEdgeCloud,downsampledSurfCloud);

}

// 根据当前欧拉角的参数更新旋转四元数和平移矩阵
void OdomEstimationClass::updatePose(){
    // 围绕x y z轴旋转的角度，对应pitch yaw roll
    double rx = paramEuler[0], ry = paramEuler[1], rz = paramEuler[2];
    // x y z方向平移的距离
    double tx = paramEuler[3], ty = paramEuler[4], tz = paramEuler[5];
    // Z * Y * X
    Eigen::Matrix3d R = (
        Eigen::AngleAxisd(rz, Eigen::Vector3d::UnitZ()) * 
        Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()) * 
        Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX())
        ).toRotationMatrix();
    Eigen::Vector3d T;
    T <<    tx, ty, tz;

    q_w_curr = Eigen::Quaterniond(R);
    t_w_curr = T;
}

// 通过当前变换矩阵，将pi变换到po
void OdomEstimationClass::pointAssociateToMap(pcl::PointXYZI const *const pi, pcl::PointXYZI *const po)
{
    Eigen::Isometry3d isoMat;
    Eigen::Vector3d point_curr(pi->x, pi->y, pi->z);
    Eigen::Vector3d point_w = q_w_curr * point_curr + t_w_curr;
    po->x = point_w.x();
    po->y = point_w.y();
    po->z = point_w.z();
    po->intensity = pi->intensity;
    //po->intensity = 1.0;
}

void OdomEstimationClass::downSamplingToMap(const pcl::PointCloud<pcl::PointXYZI>::Ptr& edge_pc_in, pcl::PointCloud<pcl::PointXYZI>::Ptr& edge_pc_out, const pcl::PointCloud<pcl::PointXYZI>::Ptr& surf_pc_in, pcl::PointCloud<pcl::PointXYZI>::Ptr& surf_pc_out){
    downSizeFilterEdge.setInputCloud(edge_pc_in);
    downSizeFilterEdge.filter(*edge_pc_out);
    downSizeFilterSurf.setInputCloud(surf_pc_in);
    downSizeFilterSurf.filter(*surf_pc_out);    
}

// 输入：特征点点云，特征点点云地图，优化问题，loss函数
void OdomEstimationClass::addEdgeCostFactor(const pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_in, const pcl::PointCloud<pcl::PointXYZI>::Ptr& map_in, ceres::Problem& problem, ceres::LossFunction *loss_function){
    int corner_num=0;
    // 遍历输入的特征点
    for (int i = 0; i < (int)pc_in->points.size(); i++)
    { 
        // 经过当前待优化的变换矩阵变换得到的点坐标
        pcl::PointXYZI point_temp;
        pointAssociateToMap(&(pc_in->points[i]), &point_temp);

        // k临近查找到的点的index
        #ifdef NANOFLANN
        std::vector<size_t> pointSearchInd;
        #else
        std::vector<int> pointSearchInd;
        #endif
        // k临近查找到的点的距离
        std::vector<float> pointSearchSqDis;
        // 搜索point_temp周围的5个点
        kdtreeEdgeMap->nearestKSearch(point_temp, 5, pointSearchInd, pointSearchSqDis); 
        // 如果最远点距离也比较近
        if (pointSearchSqDis[4] < 1.0)
        {
            // 依次储存5个最近点
            std::vector<Eigen::Vector3d> nearCorners;
            // 储存5个点的重心坐标
            Eigen::Vector3d center(0, 0, 0);
            // 取得5个点的重心
            for (int j = 0; j < 5; j++)
            {
                Eigen::Vector3d tmp(map_in->points[pointSearchInd[j]].x,
                                    map_in->points[pointSearchInd[j]].y,
                                    map_in->points[pointSearchInd[j]].z);
                center = center + tmp;
                nearCorners.push_back(tmp);
            }
            center = center / 5.0;

            // 5个点协方差矩阵
            // https://njuferret.github.io/2019/07/28/2019-07-28_geometric-interpretation-covariance-matrix/
            Eigen::Matrix3d covMat = Eigen::Matrix3d::Zero();
            for (int j = 0; j < 5; j++)
            {
                Eigen::Matrix<double, 3, 1> tmpZeroMean = nearCorners[j] - center;
                covMat = covMat + tmpZeroMean * tmpZeroMean.transpose();
            }

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(covMat);

            // 5个点构成的线的方向
            Eigen::Vector3d unit_direction = saes.eigenvectors().col(2);
            // 未经变换的当前点，和point_temp不同
            Eigen::Vector3d curr_point(pc_in->points[i].x, pc_in->points[i].y, pc_in->points[i].z);
            // 如果这5个点确实近似构成一条直线
            if (saes.eigenvalues()[2] > 3 * saes.eigenvalues()[1])
            { 
                // 5个点的中心作为在直线上的点
                Eigen::Vector3d point_on_line = center;
                Eigen::Vector3d point_a, point_b;
                // 在直线上找2个点
                point_a = 0.1 * unit_direction + point_on_line;
                point_b = -0.1 * unit_direction + point_on_line;

                // 添加一个边的误差项
                ceres::CostFunction *cost_function = new EdgeAnalyticCostFunction(curr_point, point_a, point_b);  
                problem.AddResidualBlock(cost_function, loss_function, paramEuler);
                corner_num++;   
            }                           
        }
    }
    if(corner_num<20){
        printf("not enough correct points");
    }

}

void OdomEstimationClass::addSurfCostFactor(const pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_in, const pcl::PointCloud<pcl::PointXYZI>::Ptr& map_in, ceres::Problem& problem, ceres::LossFunction *loss_function){
    int surf_num=0;
    for (int i = 0; i < (int)pc_in->points.size(); i++)
    {
        pcl::PointXYZI point_temp;
        pointAssociateToMap(&(pc_in->points[i]), &point_temp);
        #ifdef NANOFLANN
        std::vector<size_t> pointSearchInd;
        #else
        std::vector<int> pointSearchInd;
        #endif
        std::vector<float> pointSearchSqDis;
        // 在经过变换之后的面特征地图中查找最近的5个点
        kdtreeSurfMap->nearestKSearch(point_temp, 5, pointSearchInd, pointSearchSqDis);

        // 每一行储存一个点的坐标
        Eigen::Matrix<double, 5, 3> matA0;
        Eigen::Matrix<double, 5, 1> matB0 = -1 * Eigen::Matrix<double, 5, 1>::Ones();
        if (pointSearchSqDis[4] < 1.0)
        {
            // 向matA0中添加点
            for (int j = 0; j < 5; j++)
            {
                matA0(j, 0) = map_in->points[pointSearchInd[j]].x;
                matA0(j, 1) = map_in->points[pointSearchInd[j]].y;
                matA0(j, 2) = map_in->points[pointSearchInd[j]].z;
            }
            // find the norm of plane，解方程matA0 . norm = matB0， 此时norm和5个点任意一个点的内积接近-1，norm是平面法向量
            Eigen::Vector3d norm = matA0.colPivHouseholderQr().solve(matB0);
            // norm的长度倒数
            double negative_OA_dot_norm = 1 / norm.norm();
            // 单位化
            norm.normalize();

            bool planeValid = true;
            for (int j = 0; j < 5; j++)
            {
                // if OX * n > 0.2, then plane is not fit well，平面法向量和平面点的内积，如果严格垂直的话内积为0
                if (fabs(norm(0) * map_in->points[pointSearchInd[j]].x +
                         norm(1) * map_in->points[pointSearchInd[j]].y +
                         norm(2) * map_in->points[pointSearchInd[j]].z + negative_OA_dot_norm) > 0.2)
                {
                    planeValid = false;
                    break;
                }
            }
            Eigen::Vector3d curr_point(pc_in->points[i].x, pc_in->points[i].y, pc_in->points[i].z);
            // 如果平面存在
            if (planeValid)
            {
                // 添加误差项
                ceres::CostFunction *cost_function = new SurfNormAnalyticCostFunction(curr_point, norm, negative_OA_dot_norm);    
                problem.AddResidualBlock(cost_function, loss_function, paramEuler);

                surf_num++;
            }
        }

    }
    if(surf_num<20){
        printf("not enough correct points");
    }

}

// 添加特征点到edge和surf地图中
void OdomEstimationClass::addPointsToMap(const pcl::PointCloud<pcl::PointXYZI>::Ptr& downsampledEdgeCloud, const pcl::PointCloud<pcl::PointXYZI>::Ptr& downsampledSurfCloud){

    // 将cloud中经过优化完成的变换矩阵的变换，然后添加到对应的地图中
    for (int i = 0; i < (int)downsampledEdgeCloud->points.size(); i++)
    {
        pcl::PointXYZI point_temp;
        pointAssociateToMap(&downsampledEdgeCloud->points[i], &point_temp);
        laserCloudCornerMap->push_back(point_temp); 
    }
    
    for (int i = 0; i < (int)downsampledSurfCloud->points.size(); i++)
    {
        pcl::PointXYZI point_temp;
        pointAssociateToMap(&downsampledSurfCloud->points[i], &point_temp);
        laserCloudSurfMap->push_back(point_temp);
    }
    
    // 更新local map大小
    double x_min = +odom.translation().x()-100;
    double y_min = +odom.translation().y()-100;
    double z_min = +odom.translation().z()-100;
    double x_max = +odom.translation().x()+100;
    double y_max = +odom.translation().y()+100;
    double z_max = +odom.translation().z()+100;
    
    //ROS_INFO("size : %f,%f,%f,%f,%f,%f", x_min, y_min, z_min,x_max, y_max, z_max);
    cropBoxFilter.setMin(Eigen::Vector4f(x_min, y_min, z_min, 1.0));
    cropBoxFilter.setMax(Eigen::Vector4f(x_max, y_max, z_max, 1.0));
    cropBoxFilter.setNegative(false);    

    // 重新选取特征点云地图，以防特征点云地图越来越大
    pcl::PointCloud<pcl::PointXYZI>::Ptr tmpCorner(new pcl::PointCloud<pcl::PointXYZI>());
    pcl::PointCloud<pcl::PointXYZI>::Ptr tmpSurf(new pcl::PointCloud<pcl::PointXYZI>());
    cropBoxFilter.setInputCloud(laserCloudSurfMap);
    cropBoxFilter.filter(*tmpSurf);
    cropBoxFilter.setInputCloud(laserCloudCornerMap);
    cropBoxFilter.filter(*tmpCorner);

    downSizeFilterSurf.setInputCloud(tmpSurf);
    downSizeFilterSurf.filter(*laserCloudSurfMap);
    downSizeFilterEdge.setInputCloud(tmpCorner);
    downSizeFilterEdge.filter(*laserCloudCornerMap);

}

void OdomEstimationClass::getMap(pcl::PointCloud<pcl::PointXYZI>::Ptr& laserCloudMap){
    
    *laserCloudMap += *laserCloudSurfMap;
    *laserCloudMap += *laserCloudCornerMap;
}

OdomEstimationClass::OdomEstimationClass(){

}
