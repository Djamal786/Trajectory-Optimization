function out1 = T_full_3(in1,in2)
%T_full_3
%    OUT1 = T_full_3(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    29-Apr-2023 17:10:28

d1 = in2(1,:);
q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = sin(q1);
t6 = sin(q2);
t7 = sin(q3);
out1 = reshape([-t5.*t7+t2.*t3.*t4,t2.*t7+t3.*t4.*t5,t4.*t6,0.0,t2.*t6,t5.*t6,-t3,0.0,-t4.*t5-t2.*t3.*t7,t2.*t4-t3.*t5.*t7,-t6.*t7,0.0,-d1.*t2.*t6,-d1.*t5.*t6,d1.*t3,1.0],[4,4]);
