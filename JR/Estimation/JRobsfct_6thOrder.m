function r = JRobsfct_6thOrder(x)

global dq

C = [zeros(1,dq) 0 0 1 0 -1 0];           % observation function

r = C*x;
