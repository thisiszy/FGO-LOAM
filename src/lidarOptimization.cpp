// Author of FLOAM: Wang Han 
// Email wh200720041@gmail.com
// Homepage https://wanghan.pro

#include "lidarOptimization.h"

// 参数：需要优化的点的坐标，直线上点a，直线上点b
EdgeAnalyticCostFunction::EdgeAnalyticCostFunction(Eigen::Vector3d curr_point_, Eigen::Vector3d last_point_a_, Eigen::Vector3d last_point_b_)
        : curr_point(curr_point_), last_point_a(last_point_a_), last_point_b(last_point_b_){

}

bool EdgeAnalyticCostFunction::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    // 本次优化前的旋转矩阵
    Eigen::Quaterniond q_w_curr;
    // 本次优化前的平移矩阵
    Eigen::Vector3d t_w_curr;
    getPose(parameters, q_w_curr, t_w_curr);

    // 通过以上位姿变换参数变换得到的点的坐标
    Eigen::Vector3d point_w = q_w_curr * curr_point + t_w_curr; 

    // 使用公式residuals[0]=|(lp-a)x(lp-b)|/|a-b|（lp-a和lp-b两个向量构成的平行四边形面积/平行四边形对角边==点lp到a-b的距离）
    Eigen::Vector3d nu = (point_w - last_point_a).cross(point_w - last_point_b);
    Eigen::Vector3d de = last_point_b - last_point_a;
    residuals[0] = nu.norm()/de.norm();
    /*
    d(R)/d(rx) = 
    {
        {0, Cos[x] Cos[z] Sin[y] + Sin[x] Sin[z],   -Cos[z] Sin[x] Sin[y] + Cos[x] Sin[z]   }, 
        {0, -Cos[z] Sin[x] + Cos[x] Sin[y] Sin[z],  -Cos[x] Cos[z] - Sin[x] Sin[y] Sin[z]   }, 
        {0, Cos[x] Cos[y],                          -Cos[y] Sin[x]                          }
    }
    d(R)/d(ry) = 
    {
        {-Cos[z] Sin[y],    Cos[y] Cos[z] Sin[x],   Cos[x] Cos[y] Cos[z]}, 
        {-Sin[y] Sin[z],    Cos[y] Sin[x] Sin[z],   Cos[x] Cos[y] Sin[z]}, 
        {-Cos[y],           -Sin[x] Sin[y],         -Cos[x] Sin[y]      }
    }
    d(R)/d(rz) = 
    {
        {-Cos[y] Sin[z],    -Cos[x] Cos[z] - Sin[x] Sin[y] Sin[z],  Cos[z] Sin[x] - Cos[x] Sin[y] Sin[z]}, 
        {Cos[y] Cos[z],     Cos[z] Sin[x] Sin[y] - Cos[x] Sin[z],   Cos[x] Cos[z] Sin[y] + Sin[x] Sin[z]}, 
        {0,                 0,                                      0                                   }
    }
    */
    if(jacobians != NULL && jacobians[0] != NULL)
    {
        Eigen::Matrix3d dR_drx;
        Eigen::Matrix3d dR_dry;
        Eigen::Matrix3d dR_drz;
        Eigen::Vector3d dRp_drx;
        Eigen::Vector3d dRp_dry;
        Eigen::Vector3d dRp_drz;
        double rx = parameters[0][0], ry = parameters[0][1], rz = parameters[0][2];
        double tx = parameters[0][3], ty = parameters[0][4], tz = parameters[0][5];
        dR_drx<<0.0, cos(rx) * cos(rz) * sin(ry) + sin(rx) * sin(rz),   -cos(rz) * sin(rx) * sin(ry) + cos(rx) * sin(rz)    , 
                0.0, -cos(rz) * sin(rx) + cos(rx) * sin(ry) * sin(rz),  -cos(rx) * cos(rz) - sin(rx) * sin(ry) * sin(rz)    , 
                0.0, cos(rx) * cos(ry),                                 -cos(ry) * sin(rx)                                  ;
        dR_dry<<-cos(rz) * sin(ry), cos(ry) * cos(rz) * sin(rx),    cos(rx) * cos(ry) * cos(rz)                                     , 
                -sin(ry) * sin(rz), cos(ry) * sin(rx) * sin(rz),    cos(rx) * cos(ry) * sin(rz)                                     , 
                -cos(ry),           -sin(rx) * sin(ry),             -cos(rx) * sin(ry)                                              ;   
        dR_drz<<-cos(ry) * sin(rz),     -cos(rx) * cos(rz) - sin(rx) * sin(ry) * sin(rz), cos(rz) * sin(rx) - cos(rx) * sin(ry) * sin(rz)   , 
                cos(ry) * cos(rz),      cos(rz) * sin(rx) * sin(ry) - cos(rx) * sin(rz),  cos(rx) * cos(rz) * sin(ry) + sin(rx) * sin(rz)   , 
                0,                      0,                                                0                                                 ;                            
        dRp_drx = dR_drx * curr_point;
        dRp_dry = dR_dry * curr_point;
        dRp_drz = dR_drz * curr_point;
        Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor> > dL_dr(jacobians[0]);
        Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor> > dL_dt(jacobians[0]+3);

        dL_dt = (nu.cross(de) / nu.norm() / de.norm()).transpose();
        // dL_dt[0] = nu.dot(Eigen::Vector3d::UnitX().cross(de)) / nu.norm();
        // dL_dt[1] = nu.dot(Eigen::Vector3d::UnitY().cross(de)) / nu.norm();
        // dL_dt[2] = nu.dot(Eigen::Vector3d::UnitZ().cross(de)) / nu.norm();
        dL_dr[0] = dL_dt.dot(dRp_drx);
        dL_dr[1] = dL_dt.dot(dRp_dry);
        dL_dr[2] = dL_dt.dot(dRp_drz);
        
        // dL_dt[2] = 0.0;
        // dL_dr[0] = 0.0;
        // dL_dr[1] = 0.0;
    }  

    return true;
 
}   


SurfNormAnalyticCostFunction::SurfNormAnalyticCostFunction(Eigen::Vector3d curr_point_, Eigen::Vector3d plane_unit_norm_, double negative_OA_dot_norm_, double weight_) 
                                                        : curr_point(curr_point_), plane_unit_norm(plane_unit_norm_), negative_OA_dot_norm(negative_OA_dot_norm_), weight(weight_){

}

bool SurfNormAnalyticCostFunction::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    // 本次优化前的旋转矩阵
    Eigen::Quaterniond q_w_curr;
    // 本次优化前的平移矩阵
    Eigen::Vector3d t_w_curr;
    getPose(parameters, q_w_curr, t_w_curr);
    Eigen::Vector3d point_w = q_w_curr * curr_point + t_w_curr;
    // residuals[0]=平面法向量 . w的向量 + 法向量的反向，这个值越小越好
    residuals[0] = weight * (plane_unit_norm.dot(point_w) + negative_OA_dot_norm);

    if(jacobians != NULL && jacobians[0] != NULL)
    {
        Eigen::Matrix3d dR_drx;
        Eigen::Matrix3d dR_dry;
        Eigen::Matrix3d dR_drz;
        Eigen::Vector3d dRp_drx;
        Eigen::Vector3d dRp_dry;
        Eigen::Vector3d dRp_drz;
        double rx = parameters[0][0], ry = parameters[0][1], rz = parameters[0][2];
        double tx = parameters[0][3], ty = parameters[0][4], tz = parameters[0][5];
        dR_drx<<0.0, cos(rx) * cos(rz) * sin(ry) + sin(rx) * sin(rz),   -cos(rz) * sin(rx) * sin(ry) + cos(rx) * sin(rz)    , 
                0.0, -cos(rz) * sin(rx) + cos(rx) * sin(ry) * sin(rz),  -cos(rx) * cos(rz) - sin(rx) * sin(ry) * sin(rz)    , 
                0.0, cos(rx) * cos(ry),                                 -cos(ry) * sin(rx)                                  ;
        dR_dry<<-cos(rz) * sin(ry), cos(ry) * cos(rz) * sin(rx),    cos(rx) * cos(ry) * cos(rz)                                     , 
                -sin(ry) * sin(rz), cos(ry) * sin(rx) * sin(rz),    cos(rx) * cos(ry) * sin(rz)                                     , 
                -cos(ry),           -sin(rx) * sin(ry),             -cos(rx) * sin(ry)                                              ;   
        dR_drz<<-cos(ry) * sin(rz),     -cos(rx) * cos(rz) - sin(rx) * sin(ry) * sin(rz), cos(rz) * sin(rx) - cos(rx) * sin(ry) * sin(rz)   , 
                cos(ry) * cos(rz),      cos(rz) * sin(rx) * sin(ry) - cos(rx) * sin(rz),  cos(rx) * cos(rz) * sin(ry) + sin(rx) * sin(rz)   , 
                0,                      0,                                                0                                                 ;                            
        dRp_drx = dR_drx * curr_point;
        dRp_dry = dR_dry * curr_point;
        dRp_drz = dR_drz * curr_point;
        Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor> > dL_dr(jacobians[0]);
        Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor> > dL_dt(jacobians[0]+3);
        dL_dt = plane_unit_norm.transpose();
        dL_dr[0] = dL_dt.dot(dRp_drx);
        dL_dr[1] = dL_dt.dot(dRp_dry);
        dL_dr[2] = dL_dt.dot(dRp_drz);

        // dL_dr[2] = 0.0;
        // dL_dt[0] = 0.0;
        // dL_dt[1] = 0.0;
    }
    return true;

}   

void getPose(double const * const *paramEuler, Eigen::Quaterniond& q, Eigen::Vector3d& t){
    // 围绕x y z轴旋转的角度，对应pitch yaw roll
    double rx = paramEuler[0][0], ry = paramEuler[0][1], rz = paramEuler[0][2];
    // x y z方向平移的距离
    double tx = paramEuler[0][3], ty = paramEuler[0][4], tz = paramEuler[0][5];
    // Z * Y * X
    Eigen::Matrix3d R = (
        Eigen::AngleAxisd(rz, Eigen::Vector3d::UnitZ()) * 
        Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()) * 
        Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX())
        ).toRotationMatrix();
    Eigen::Vector3d T;
    T <<    tx, ty, tz;

    q = Eigen::Quaterniond(R);
    t = T;
}
