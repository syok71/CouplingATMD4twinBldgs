% Augmanr Mass Matrix

function M=addatmd0(n1,M1,n2,M2)

M=[M1 zeros(n1,n2);zeros(n2,n1) M2];