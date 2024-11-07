/*This file is part of lbr_iiwa_7_R800 package.
 * Copyright (C) 2019, Mario Selvaggio

 * Email id : mario.selvaggio@unina.it
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *   * Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the names of Stanford University or Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
* This code will subscriber integer values from demo_topic_publisher
*/

#include <stdio.h>
#include <math.h>
//#include <lbr_iiwa_7_kinematics_pkg/lbr_iiwa_7_robot.h>
#include "iiwa_qp/lbr_iiwa_7_robot.h"

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

//******************************************************************************
lbr_iiwa_7_robot::lbr_iiwa_7_robot(const KDL::JntArray &_q_min,
                                   const KDL::JntArray &_q_max)
{
    q_min_ = new KDL::JntArray(_q_min);
    q_max_ = new KDL::JntArray(_q_max);
}

lbr_iiwa_7_robot::lbr_iiwa_7_robot(){}
//******************************************************************************

//Forward Kinematics Matrix Te
KDL::Frame lbr_iiwa_7_robot::forwardKinematics(const KDL::JntArray &q){
    KDL::Frame A0;

    double q1 = q.data[0];
    double q2 = q.data[1];
    double q3 = q.data[2];
    double q4 = q.data[3];
    double q5 = q.data[4];
    double q6 = q.data[5];
    double q7 = q.data[6];

    double t2 = cos(q1);
    double t3 = cos(q4);
    double t4 = sin(q1);
    double t5 = sin(q3);
    double t6 = t4*t5;
    double t7 = cos(q2);
    double t8 = cos(q3);
    double t13 = t2*t7*t8;
    double t9 = t6-t13;
    double t10 = sin(q2);
    double t11 = sin(q4);
    double t12 = sin(q5);
    double t14 = t3*t9;
    double t27 = t2*t10*t11;
    double t15 = t14-t27;
    double t16 = cos(q5);
    double t17 = t4*t8;
    double t18 = t2*t5*t7;
    double t19 = t17+t18;
    double t20 = sin(q7);
    double t21 = sin(q6);
    double t22 = t9*t11;
    double t23 = t2*t3*t10;
    double t24 = t22+t23;
    double t25 = t21*t24;
    double t26 = cos(q6);
    double t28 = t15*t16;
    double t29 = t12*t19;
    double t30 = t28+t29;
    double t31 = t26*t30;
    double t32 = t25+t31;
    double t33 = cos(q7);
    double t34 = t12*t15;
    double t35 = t34-t16*t19;
    double t36 = t2*t5;
    double t37 = t4*t7*t8;
    double t38 = t36+t37;
    double t39 = t3*t38;
    double t40 = t4*t10*t11;
    double t41 = t39+t40;
    double t42 = t2*t8;
    double t48 = t4*t5*t7;
    double t43 = t42-t48;
    double t44 = t11*t38;
    double t55 = t3*t4*t10;
    double t45 = t44-t55;
    double t46 = t21*t45;
    double t47 = t16*t41;
    double t49 = t12*t43;
    double t50 = t47+t49;
    double t51 = t26*t50;
    double t52 = t46+t51;
    double t53 = t12*t41;
    double t54 = t53-t16*t43;
    double t56 = t7*t11;
    double t58 = t3*t8*t10;
    double t57 = t56-t58;
    double t59 = t12*t57;
    double t60 = t59-t5*t10*t16;
    double t61 = t16*t57;
    double t62 = t5*t10*t12;
    double t63 = t61+t62;
    double t64 = t3*t7;
    double t65 = t8*t10*t11;
    double t66 = t64+t65;
    double t67 = t21*t66;
    A0.M(0,0) = t20*t35-t32*t33;
    A0.M(0,1) = t20*t32+t33*t35;
    A0.M(0,2) = t24*t26-t21*t30;
    A0.p(0) = t2*t10*(2.0/5.0)+t9*t11*(2.0/5.0)+t24*t26*(1.11E2/1.0E3)-t21*t30*(1.11E2/1.0E3)+t2*t3*t10*(2.0/5.0);
    A0.M(1,0) = -t20*t54+t33*t52;
    A0.M(1,1) = -t20*t52-t33*t54;
    A0.M(1,2) = t21*t50-t26*t45;
    A0.p(1) = t4*t10*(2.0/5.0)-t11*t38*(2.0/5.0)+t21*t50*(1.11E2/1.0E3)-t26*t45*(1.11E2/1.0E3)+t3*t4*t10*(2.0/5.0);
    A0.M(2,0) = -t20*t60-t33*(t67-t26*t63);
    A0.M(2,1) = -t33*t60+t20*(t67-t26*t63);
    A0.M(2,2) = t21*t63+t26*t66;
    A0.p(2) = t7*(2.0/5.0)+t3*t7*(2.0/5.0)+t21*t63*(1.11E2/1.0E3)+t26*t66*(1.11E2/1.0E3)+t8*t10*t11*(2.0/5.0)+1.7E1/5.0E1;

    return A0;
}

//Jacobian matrix J
KDL::Jacobian lbr_iiwa_7_robot::jacobianMatrix(const KDL::JntArray &q){

    KDL::Jacobian A0(7);
    A0.data.setZero();

    double q1 = q.data[0];
    double q2 = q.data[1];
    double q3 = q.data[2];
    double q4 = q.data[3];
    double q5 = q.data[4];
    double q6 = q.data[5];
    double q7 = q.data[6];

    double t2 = sin(q1);
    double t3 = sin(q2);
    double t4 = cos(q4);
    double t5 = cos(q1);
    double t6 = sin(q3);
    double t7 = t5*t6;
    double t8 = cos(q2);
    double t9 = cos(q3);
    double t10 = t2*t8*t9;
    double t11 = t7+t10;
    double t12 = sin(q4);
    double t13 = sin(q6);
    double t14 = cos(q5);
    double t15 = sin(q5);
    double t16 = cos(q6);
    double t17 = t11*t12;
    double t45 = t2*t3*t4;
    double t18 = t17-t45;
    double t19 = t16*t18*(1.11E2/1.0E3);
    double t20 = t4*t11;
    double t21 = t2*t3*t12;
    double t22 = t20+t21;
    double t23 = t14*t22;
    double t24 = t5*t9;
    double t43 = t2*t6*t8;
    double t25 = t24-t43;
    double t26 = t15*t25;
    double t27 = t23+t26;
    double t28 = t11*t12*(2.0/5.0);
    double t29 = t8*(2.0/5.0);
    double t30 = t4*t8*(2.0/5.0);
    double t31 = t8*t12;
    double t44 = t3*t4*t9;
    double t32 = t31-t44;
    double t33 = t14*t32;
    double t34 = t3*t6*t15;
    double t35 = t33+t34;
    double t36 = t13*t35*(1.11E2/1.0E3);
    double t37 = t4*t8;
    double t38 = t3*t9*t12;
    double t39 = t37+t38;
    double t40 = t16*t39*(1.11E2/1.0E3);
    double t41 = t3*t9*t12*(2.0/5.0);
    double t42 = t29+t30+t36+t40+t41;
    double t46 = t13*t27*(1.11E2/1.0E3);
    double t47 = t2*t3*t4*(2.0/5.0);
    double t48 = t19+t28-t46-t47;
    double t49 = t30+t36+t40+t41;
    double t50 = t19-t46;
    double t51 = t36+t40;
    double t52 = t2*t6;
    double t54 = t5*t8*t9;
    double t53 = t52-t54;
    double t55 = t3*t5*(2.0/5.0);
    double t56 = t12*t53;
    double t57 = t3*t4*t5;
    double t58 = t56+t57;
    double t59 = t16*t58*(1.11E2/1.0E3);
    double t60 = t4*t53;
    double t71 = t3*t5*t12;
    double t61 = t60-t71;
    double t62 = t14*t61;
    double t63 = t2*t9;
    double t64 = t5*t6*t8;
    double t65 = t63+t64;
    double t66 = t15*t65;
    double t67 = t62+t66;
    double t68 = t12*t53*(2.0/5.0);
    double t69 = t3*t4*t5*(2.0/5.0);
    double t72 = t13*t67*(1.11E2/1.0E3);
    double t70 = t55+t59+t68+t69-t72;
    double t73 = t59+t68+t69-t72;
    double t74 = t15*t32;
    double t94 = t3*t6*t14;
    double t75 = t74-t94;
    double t76 = t13*t35;
    double t77 = t16*t39;
    double t78 = t76+t77;
    double t79 = t59-t72;
    double t80 = t2*t3*(2.0/5.0);
    double t81 = -t19-t28+t46+t47+t80;
    double t82 = t15*t22;
    double t92 = t14*t25;
    double t83 = t82-t92;
    double t84 = t15*t61;
    double t90 = t14*t65;
    double t85 = t84-t90;
    double t86 = t16*t18;
    double t93 = t13*t27;
    double t87 = t86-t93;
    double t88 = t16*t58;
    double t91 = t13*t67;
    double t89 = t88-t91;
    A0(0,0) = t19+t28-t2*t3*(2.0/5.0)-t13*t27*(1.11E2/1.0E3)-t2*t3*t4*(2.0/5.0);
    A0(0,1) = t5*t42;
    A0(0,2) = -t8*t81+t2*t3*t42;
    A0(0,3) = -t25*t49-t3*t6*t48;
    A0(0,4) = -t18*t49+t39*t48;
    A0(0,5) = -t50*t75-t51*t83;
    A0(0,6) = t50*t78-t51*t87;
    A0(1,0) = t70;
    A0(1,1) = t2*t42;
    A0(1,2) = t8*t70-t3*t5*t42;
    A0(1,3) = -t49*t65-t3*t6*t73;
    A0(1,4) = -t49*t58+t39*t73;
    A0(1,5) = -t51*t85-t75*t79;
    A0(1,6) = -t51*t89+t78*t79;
    A0(2,1) = -t5*t70-t2*t81;
    A0(2,2) = -t2*t3*t70+t3*t5*t81;
    A0(2,3) = t25*t73-t48*t65;
    A0(2,4) = t18*t73-t48*t58;
    A0(2,5) = -t50*t85+t79*t83;
    A0(2,6) = -t50*t89+t79*t87;
    A0(3,1) = -t2;
    A0(3,2) = t3*t5;
    A0(3,3) = t65;
    A0(3,4) = t58;
    A0(3,5) = t85;
    A0(3,6) = t89;
    A0(4,1) = t5;
    A0(4,2) = t2*t3;
    A0(4,3) = -t24+t43;
    A0(4,4) = -t17+t45;
    A0(4,5) = -t82+t92;
    A0(4,6) = -t86+t93;
    A0(5,0) = 1.0;
    A0(5,2) = t8;
    A0(5,3) = -t3*t6;
    A0(5,4) = t39;
    A0(5,5) = -t74+t94;
    A0(5,6) = t78;

    return A0;
}

// inverse kinematics
KDL::JntArray lbr_iiwa_7_robot::backwardKinematics(const KDL::Frame& T_d,
                                            const KDL::JntArray& q_i)
{
    // std::cout << "calculating inverse kinematics... " << std::endl;
    double _Tsam = 0.001;        // sampling time
    double _max_pError = 0.0001;  // max norm position error
    double _max_oError = 0.0005;  // max norm orientation error
    double _Kp = 1000;           // position gain
    double _Ko = 500;            // orientation gain
    int _n_max_iter = 100;       // maximum iteration number
    Matrix3d Kp = _Kp*Matrix3d::Identity();
    Matrix3d Ko = _Ko*Matrix3d::Identity();

    Vector6d q;
    q << q_i.data;
    KDL::JntArray q_kdl(6);
    q_kdl.data << q;

    Vector6d dq = Vector6d::Zero();
    Vector6d dp = Vector6d::Zero();
    Matrix6d W_inv = Matrix6d::Identity();

    KDL::Frame T_e = forwardKinematics(q_kdl);
    KDL::Jacobian J = jacobianMatrix(q_kdl);

    Matrix6d pinvJ = Matrix6d::Zero();
    Matrix6d I = Matrix6d::Identity();

    Matrix6d J_T = Matrix6d::Zero();
    Matrix3d L = Matrix3d::Zero();
    Matrix3d pinvL = Matrix3d::Zero();
    Vector3d vp = Vector3d::Zero();
    Vector3d vo = Vector3d::Zero();
    Vector6d v = Vector6d::Zero();

    int n_iter = 0;

    Vector3d p_d(T_d.p.data);
    Vector3d p_e(T_e.p.data);
    Vector3d ep;
    ep = p_d - p_e;

    Eigen::Matrix<double,3,3,RowMajor> Rd(T_d.M.data);
    Eigen::Matrix<double,3,3,RowMajor> Re(T_e.M.data);
    Vector3d eo;
    eo = utilities::rotationMatrixError(Rd, Re);
    double pError = ep.norm();                                     //L_2 norm of vector
    double oError = eo.norm();

    while(n_iter<=_n_max_iter && (pError>_max_pError || oError>_max_oError)){

        //L = utilities::L_matrix(Td.block(0,0,3,3), Te.block(0,0,3,3));
        //FullPivLU<Matrix3d> Core_L(L);
        //pinvL = Core_L.inverse();

        //vp = dp.block(0,0,3,1) + Kp*ep;
        vp =  Kp*ep;

        //vo = pinvL*(L.transpose()*dp.block(3,0,3,1) + Ko*eo);
        vo = Ko*eo;
        v << vp[0], vp[1], vp[2], vo[0], vo[1], vo[2];
        // std::cout << vo << std::endl;
        //J_T = J.transpose();

        // If we need a weighted pinv
        //FullPivLU<Matrix6d> Core_J(J*W_inv*J_T);
        //pinvJ = W_inv*J_T*Core_J.inverse();
        //FullPivLU<Matrix6d> Core_J(J);

        J = jacobianMatrix(q_kdl);
        const Matrix<double,6,6> j_eigen = J.data.block(0,0,6,6);
        FullPivLU<Matrix6d> lu(j_eigen);
        lu.setThreshold(1e-5);
        if (lu.rank() < 6)
        {
            std::cout << "[WARNING]: jacobian rank < 6!" << std::endl;
            break;
        }

        //pinvJ = J.data.inverse();
        pinvJ = J.data.transpose()*(J.data*J.data.transpose()+0.00001*Eigen::Matrix<double,6,6>::Identity()).inverse();
        //pinvJ = J.data.transpose()*(J.data*J.data.transpose() + 0.001*I).inverse();

        dq = pinvJ*v;
        q += dq*_Tsam;
        // std::cout << q.transpose() << std::endl;
        for(unsigned int i = 0; i < 6; i++)
        {
            if (q[i] < -M_PI)
            {
                std::cout << "[WARNING] joint " << i << " < - PI \n";
                q[i] = M_PI + (q[i] + M_PI);
            }

            if (q[i] > M_PI)
            {
                std::cout << "[WARNING] joint " << i << " > PI \n";
                q[i] = - M_PI + (q[i] - M_PI);
            }
        }

        for(unsigned int i = 0; i < 6; i++)
        {
            if (q[i] < q_min_->data[i]+0.001)
            {
                std::cout << "[WARNING]: joint " << i << " < lower limit. Value = "<< q(i,0) << "\n";
                q[i] = q_min_->data[i]+0.001;
            }
            if (q[i] > q_max_->data[i]-0.001)
            {
                std::cout << "[WARNING]: joint " << i << " > upper limit. Value = "<< q(i,0) << "\n";
                q[i] = q_max_->data[i]-0.001;
            }
        }

        q_kdl.data << q;

        T_e = forwardKinematics(q_kdl);

        p_e = Vector3d(T_e.p.data);
        Re = Eigen::Matrix<double,3,3,RowMajor>(T_e.M.data);

        ep = p_d - p_e;

        eo = utilities::rotationMatrixError(Rd, Re);

        pError = ep.norm(); //L_2 norm of vector
        oError = eo.norm();
        n_iter++;

//        std::cout << "q : " << q.transpose() << std::endl;
//        std::cout << "q_kdl : " << q_kdl.data.transpose() << std::endl;
//        std::cout << "T_d : " << T_d << std::endl;
//        std::cout << "T_e : " << T_e << std::endl;
//        std::cout << "J : " << J.data << std::endl;
//        std::cout << "R_d : " << Rd << std::endl;
//        std::cout << "R_e : " << Re << std::endl;
//        std::cout << "eo : " << eo << std::endl;
//        std::cout << "pError: " << pError << std::endl;
//        std::cout << "oError : " << oError << std::endl;
//        int c=getchar();

    }

//    std::cout << "q : " << q_kdl.data.transpose() << std::endl;
//    std::cout << "Re : " << Re << std::endl;
//    std::cout << "Rd : " << Rd << std::endl;
//    std::cout << "# iter: " << n_iter << std::endl;
//    cout << "position error: " << pError << std::endl
//         << "orientation Error: " << oError << std::endl;

    if (n_iter >= _n_max_iter)
    {
        std::cout << "[WARNING]: Ik maximum iterations reached." << std::endl;
        cout << "position error: " << pError << std::endl
             << "orientation Error: " << oError << std::endl;
    }
    q_kdl.data << q;
    return q_kdl;
}

void lbr_iiwa_7_robot::getJointLimits(KDL::JntArray& q_max,
                                      KDL::JntArray& q_min)
{
    q_max.data = this->q_max_->data;
    q_min.data = this->q_min_->data;
}


