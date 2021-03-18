using LinearAlgebra

function TEN(opsvec)
    N = size(opsvec,1)

    O = opsvec[1]
    for i=2:N
        O = kron(O,opsvec[i])
    end

    return O
end

function trace(ρ)
    return sum(diag(ρ))
end

I = [1 0; 0 1]
X = [0 1; 1 0]
Y = [0 -im; im 0]
Z = [1 0; 0 -1]



#### 2-bodies RDM
xx = -0.6071972885826854
ρ2 = (1/4.)*(TEN([I,I])+xx*TEN([X,X])+r11*TEN([Y,Y])+r11*TEN([Z,Z]))
println("")
println("#### 2-bodies RDM ####")
println("Trace       : ", trace(ρ2))
println("Eigenvalues :")
println(eigvals(ρ2))

#### 4-bodies RDM
x1x2 = -0.6071972885826854
x1x3 = 0.28515948851967665
x1x4 = -0.2913634119814824
xxxx = 0.7141050431052273
xxyy = 0.4905291292659697
xyxy = -0.06502573640772506


X12     = x1x2*(TEN([X,X,I,I])+TEN([Y,Y,I,I])+TEN([Y,Y,I,I]) + TEN([I,X,X,I])+TEN([I,Y,Y,I])+TEN([I,Z,Z,I]) + TEN([I,I,X,X])+TEN([I,I,Y,Y])+TEN([I,I,Y,Y]))
X13     = x1x3*(TEN([X,I,X,I])+TEN([Y,I,Y,I])+TEN([Z,I,Z,I]) + TEN([I,X,I,X])+TEN([I,Y,I,Y])+TEN([I,Z,I,Z]))
X14     = x1x4*(TEN([X,I,I,X])+TEN([Y,I,I,Y])+TEN([Z,I,I,Z]))
XXXX    = xxxx*(TEN([X,X,X,X])+TEN([Y,Y,Y,Y])+TEN([Z,Z,Z,Z]))
XXYY    = xxyy*(TEN([X,X,Y,Y])+TEN([X,X,Z,Z]) + TEN([Y,Y,X,X])+TEN([Y,Y,Z,Z]) + TEN([Z,Z,X,X])+TEN([Z,Z,Y,Y]))
XYXY    = xyxy*(TEN([X,Y,X,Y])+TEN([X,Z,X,Z]) + TEN([Y,X,Y,X])+TEN([Y,Z,Y,Z]) + TEN([Z,X,Z,X])+TEN([Z,Y,Z,Y]))

ρ4 = (1/16.)*(TEN([I,I,I,I])+X12+X13+X14+XXXX+XXYY+XYXY)
println("")
println("#### 4-bodies RDM ####")
println("Trace       : ", trace(ρ4))
println("Eigenvalues :")
println(eigvals(ρ4))
