#include "lec5_hw/visualizer.hpp"
#include "lec5_hw/trajectory.hpp"

#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>

#include <cmath>
#include <iostream>
#include <vector>

//using namespace std;
//using namespace Eigen;

struct Config
{
    std::string targetTopic;
    double clickHeight;
    std::vector<double> initialVel;
    std::vector<double> initialAcc;
    std::vector<double> terminalVel;
    std::vector<double> terminalAcc;
    double allocationSpeed;
    double allocationAcc;
    int maxPieceNum;

    Config(const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("ClickHeight", clickHeight);
        nh_priv.getParam("InitialVel", initialVel);
        nh_priv.getParam("InitialAcc", initialAcc);
        nh_priv.getParam("TerminalVel", terminalVel);
        nh_priv.getParam("TerminalAcc", terminalAcc);
        nh_priv.getParam("AllocationSpeed", allocationSpeed);
        nh_priv.getParam("AllocationAcc", allocationAcc);
        nh_priv.getParam("MaxPieceNum", maxPieceNum);
    }
};

double timeTrapzVel(const double dist,
                    const double vel,
                    const double acc)
{
    const double t = vel / acc;
    const double d = 0.5 * acc * t * t;

    if (dist < d + d)
    {
        return 2.0 * sqrt(dist / acc);
    }
    else
    {
        return 2.0 * t + (dist - 2.0 * d) / vel;
    }
}
/*
 * We use BIVP method:
 * general 5th order polynomial independently; minimum jerk
 * param:
 * pieceNum: the number of part of trajectory. num_piece = waypoint - 1;
 * intermediatePositions: the intermediate positions
 * timeAllocationVector: the vector of the allocation time of every piece of trajectory
 */
void minimumJerkTrajGen(
    // Inputs:
    const int pieceNum,
    const Eigen::Vector3d &initialPos,
    const Eigen::Vector3d &initialVel,
    const Eigen::Vector3d &initialAcc,
    const Eigen::Vector3d &terminalPos,
    const Eigen::Vector3d &terminalVel,
    const Eigen::Vector3d &terminalAcc,
    const Eigen::Matrix3Xd &intermediatePositions,
    const Eigen::VectorXd &timeAllocationVector,
    // Outputs:
    Eigen::MatrixX3d &coefficientMatrix) {
    // coefficientMatrix is a matrix with 6*piece num rows and 3 columes
    // As for a polynomial c0+c1*t+c2*t^2+c3*t^3+c4*t^4+c5*t^5,
    // each 6*3 sub-block of coefficientMatrix is
    // --              --
    // | c0_x c0_y c0_z |
    // | c1_x c1_y c1_z |
    // | c2_x c2_y c2_z |
    // | c3_x c3_y c3_z |
    // | c4_x c4_y c4_z |
    // | c5_x c5_y c5_z |
    // --              --
    // Please computed coefficientMatrix of the minimum-jerk trajectory
    // in this function

    // ------------------------ Put your solution below ------------------------
    // M*c = b -> (dimension * dimension)  *  (dimension * 3) = dimension * 3
    int dimen = timeAllocationVector.size()*6;
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(dimen,dimen);
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(dimen,3);
    Eigen::MatrixX3d coefficientMatrix_bak = Eigen::MatrixXd::Zero(dimen,3);;/////////////////////////
    //First,We should calculate the initial: F0, the initial part of b;
    //Then, we should calculate Ei and Fi
    // F0 of M : R^(s*2s)
    //M.block(startRow, startCol,blockRows,blockCols)
    M.block(0,0,3,6) << 1,0,0,0,0,0,\
                                                        0,1,0,0,0,0,\
                                                        0,0,2,0,0,0;
    // b: initial part 3*3
    b.block(0,0,3,3) << initialPos(0),initialPos(1),initialPos(2),\
                                                    initialVel(0),initialVel(1),initialVel(2),\
                                                    initialAcc(0),initialAcc(1),initialAcc(2);
    // Ei and Fi : R^(2s*2s) = R^(6*6)
    for(int i = 0; i < dimen/6 - 1;i++){
        int index = i*6;
        double t = timeAllocationVector(i);
        // Fi = (0,-B(0),...,-B^(di-1)(0))^T
        // B(x) is the basis Ci
        M.block(index + 3, index+6,6,6)<<0,0,0,0,0,0,\
                                    -1,0,0,0,0,0,\
                                    0,-1,0,0,0,0,\
                                    0,0,-2,0,0,0,\
                                    0,0,0,-6,0,0,\
                                    0,0,0,0,-24,0;
        b.block(index+3,0,1,3) = intermediatePositions.col(i).transpose();
        // Ei = (B(Ti),...,B^(di-1)(Ti),
        //       B(Ti),...,B^(di~-1)(Ti))^T
        //       di~ = 2s - di; di = 1
        M.block(index+3,index,6,6)<< 1, t,pow(t,2),pow(t,3),pow(t,4),pow(t,5),\
                                     1, t,pow(t,2),pow(t,3),pow(t,4),pow(t,5),\
                                     0,1,2*t,3*pow(t,2),4*pow(t,3),5*pow(t,4),\
                                     0,0,2,6*t,12*pow(t,2),20*pow(t,3),\
                                     0,0,0,6,24*t,60*pow(t,2),\
                                     0,0,0,0,24,120*t;
    }
    //EM= (B(TM),...,B^(S-1)(TM))^T
    double t = timeAllocationVector(dimen/6-1);
    M.block(dimen-3,dimen-6,3,6)<<1,t,pow(t,2),pow(t,3),pow(t,4),pow(t,5),\
                                                                0,1,2*t,3*pow(t,2),4*pow(t,3),5*pow(t,4),
                                                                0,0,2,6*t,12*pow(t,2),20*pow(t,3);
    b.block(dimen-3,0,3,3)<<terminalPos(0),terminalPos(1),terminalPos(2),\
                                                        terminalVel(0),terminalVel(1),terminalVel(2),\
                                                        terminalAcc(0),terminalAcc(1),terminalAcc(2);
    std::cout<<"M is "<<M<<std::endl;
    std::cout<<"b is "<<b<<std::endl;

    // cofficientMatrix:
    // Method: QR or lu
    // reference: https://blog.csdn.net/weixin_38009585/article/details/82355889
    clock_t time_stt = clock();
    for (int i=0;i<3;i++)
    {
        coefficientMatrix_bak.col(i)=M.colPivHouseholderQr().solve(b.col(i));
        //coefficientMatrix.col(i)=M.lu().solve(b.col(i));
    }
    std::cout<<"c (QR) is "<<coefficientMatrix_bak<<std::endl;
    std::cout<<"QR use time: " << 1000*(clock() - time_stt)/(double)CLOCKS_PER_SEC << "MS" << std::endl;

    clock_t time_stt2 = clock();
    for (int i=0;i<3;i++)
    {
        //coefficientMatrix_bak.col(i)=M.colPivHouseholderQr().solve(b.col(i));
        coefficientMatrix.col(i)=M.lu().solve(b.col(i));
    }
    std::cout<<"c (lu) is "<<coefficientMatrix<<std::endl;
    std::cout<<"LU use time: " << 1000*(clock() - time_stt2)/(double)CLOCKS_PER_SEC << "MS" << std::endl;
    // ------------------------ Put your solution above ------------------------


}

class ClickGen
{
private:
    Config config;

    ros::NodeHandle nh;
    ros::Subscriber targetSub;

    Visualizer visualizer;

    Eigen::Matrix3Xd positions;
    Eigen::VectorXd times;
    int positionNum;
    Trajectory<5> traj;

public:
    ClickGen(const Config &conf,
             ros::NodeHandle &nh_)
        : config(conf),
          nh(nh_),
          visualizer(nh),
          positions(3, config.maxPieceNum + 1),
          times(config.maxPieceNum),
          positionNum(0)
    {
        targetSub = nh.subscribe(config.targetTopic, 1,
                                 &ClickGen::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
    }

    void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {
        if (positionNum > config.maxPieceNum)
        {
            positionNum = 0;
            traj.clear();
        }

        positions(0, positionNum) = msg->pose.position.x;
        positions(1, positionNum) = msg->pose.position.y;
        positions(2, positionNum) = std::fabs(msg->pose.orientation.z) * config.clickHeight;

        if (positionNum > 0)
        {
            const double dist = (positions.col(positionNum) - positions.col(positionNum - 1)).norm();
            times(positionNum - 1) = timeTrapzVel(dist, config.allocationSpeed, config.allocationAcc);
        }

        ++positionNum;

        if (positionNum > 1)
        {
            const int pieceNum = positionNum - 1;
            const Eigen::Vector3d initialPos = positions.col(0);
            const Eigen::Vector3d initialVel(config.initialVel[0], config.initialVel[1], config.initialVel[2]);
            const Eigen::Vector3d initialAcc(config.initialAcc[0], config.initialAcc[1], config.initialAcc[2]);
            const Eigen::Vector3d terminalPos = positions.col(pieceNum);
            const Eigen::Vector3d terminalVel(config.terminalVel[0], config.terminalVel[1], config.terminalVel[2]);
            const Eigen::Vector3d terminalAcc(config.terminalAcc[0], config.terminalAcc[1], config.terminalAcc[2]);
            const Eigen::Matrix3Xd intermediatePositions = positions.middleCols(1, pieceNum - 1);
            const Eigen::VectorXd timeAllocationVector = times.head(pieceNum);

            Eigen::MatrixX3d coefficientMatrix = Eigen::MatrixXd::Zero(6 * pieceNum, 3);

            minimumJerkTrajGen(pieceNum,
                               initialPos, initialVel, initialAcc,
                               terminalPos, terminalVel, terminalAcc,
                               intermediatePositions,
                               timeAllocationVector,
                               coefficientMatrix);

            traj.clear();
            traj.reserve(pieceNum);
            for (int i = 0; i < pieceNum; i++)
            {
                traj.emplace_back(timeAllocationVector(i),
                                  coefficientMatrix.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
            }
        }

        visualizer.visualize(traj, positions.leftCols(positionNum));

        return;
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "click_gen_node");
    ros::NodeHandle nh_;
    ClickGen clickGen(Config(ros::NodeHandle("~")), nh_);
    ros::spin();
    return 0;
}
