OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
cx q[8],q[9];
cx q[7],q[6];
h q[7];
t q[8];
t q[6];
t q[7];
cx q[6],q[8];
cx q[7],q[6];
cx q[8],q[7];
tdg q[6];
cx q[8],q[6];
tdg q[8];
tdg q[6];
t q[7];
cx q[7],q[6];
cx q[8],q[7];
cx q[6],q[8];
h q[7];
cx q[7],q[5];
x q[7];
h q[7];
t q[9];
t q[5];
t q[7];
cx q[5],q[9];
cx q[7],q[5];
cx q[9],q[7];
tdg q[5];
cx q[9],q[5];
tdg q[9];
tdg q[5];
t q[7];
cx q[7],q[5];
cx q[9],q[7];
cx q[5],q[9];
h q[7];