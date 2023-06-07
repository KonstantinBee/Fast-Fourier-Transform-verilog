module FFT #(parameter Bus = 8, Bit_depth = 4) (input [Bus*Bit_depth-1:0] A, output [2*Bus*Bit_depth-1:0] B);

wire [Bus*Bit_depth-1:0] S;
wire [2*Bus*Bit_depth-1:0] D;
wire [2*Bus*Bit_depth-1:0]SD[$clog2(Bus)-1:0];

BRA #( .Bus(Bus), .Bit_depth(Bit_depth)) BRA ( .A(A), .B(S));
First_cascade # ( .Bit_depth(Bit_depth)) First_cascade [((Bus/2)-1):0] ( .X(S), .Y(D));

assign SD[0] = D;  
assign  B = SD[$clog2(Bus)-1];

genvar j;

generate
    for ( j=0 ; j < $clog2(Bus)-1; j = j+1 ) begin: FG 
    wire [2*Bus*Bit_depth-1:0]SB;
	FFT_block #( .Bus(Bus), .Bit_depth(Bit_depth+j), .Jn(j)) FFT_block_A ( .A(SD[j]), .B(SB));
    R #( .N(2*Bus*(Bit_depth+j))) R0 ( .A(SB), .B(SD[j+1]));
    end
endgenerate

endmodule

module FFT_block #(parameter Bus = 8, Bit_depth = 4, Jn = 0) (input [2*Bus*Bit_depth-1:0] A, output [2*Bus*Bit_depth-1:0] B);

wire [2*Bus*Bit_depth-1:0] W;
wire [2*Bus*Bit_depth-1:0] U;
wire [2*Bus*Bit_depth-1:0] G;
wire [2*Bus*Bit_depth-1:0] H;

Permutation_conductor  #( .Bus(Bus), .Bit_depth(Bit_depth*2), .Jn(Jn)) Permutation_conductor_A  ( .A(A), .B(W));

genvar i;

generate
    for ( i = 0; i < (Bus/2); i = i+1) 
		begin: FC 		
        Butterfly # ( .Jn(Jn), .In(i), .Bit_depth(Bit_depth)) Butterfly_A ( .YZ(W[Bit_depth*4*(i + 1) - 1 : Bit_depth*4*i]), .XAXB(U[Bit_depth*4*(i + 1) - 1 : Bit_depth*4*i]));             
        end
endgenerate

Permutation_conductor_two  #( .Bus(Bus), .Bit_depth(Bit_depth*2), .Jn(Jn+1)) Permutation_conductor_B ( .A(U), .B(B));

endmodule

module BRA #(parameter Bus = 8, Bit_depth = 2) (input [Bus*Bit_depth-1:0] A, output reg [Bus*Bit_depth-1:0] B);

//Bit reverse addressing

function [$clog2(Bus)-1:0] permutation (input [$clog2(Bus)-1:0] I);

    integer i;

    for (i = 0; i < $clog2(Bus); i = i+1) begin

       permutation [i] =  I [($clog2(Bus)-1)-i];
        
    end

endfunction

integer i;
integer j;

always @(*) begin
    for (i = 0 ; i < Bus; i = i +1) begin
        for (j = 0; j < Bit_depth ; j = j +1) begin
            B[i*Bit_depth+j] = A[permutation(i)*Bit_depth+j];
        end
    end  
end
endmodule

module First_cascade # (parameter Bit_depth = 8) (input [Bit_depth*2-1:0] X, output [4*Bit_depth-1:0] Y);

// 8/24
wire [Bit_depth-1:0] XA;
wire [Bit_depth-1:0] XB;
wire [Bit_depth-1:0] XAM;
wire [Bit_depth-1:0] XBD;
wire [Bit_depth-1:0] XBM;
reg [Bit_depth-1:0] XBDA;

KN #( .Bit_depth(2*Bit_depth)) KN_1 ( .A(X), .B(XA), .C(XB));

assign XAM = XA+XB;

integer i;

always @(*) 
    begin          
        for (i = 0; i <= Bit_depth-1; i=i+1)
        XBDA[i] <= ~XB[i];                   
    end

assign XBD = XBDA +1'd1;
assign XBM = XA + XBD; 

wire [Bit_depth-1:0] A = 0;

assign Y = {A, XBM, A, XAM};

endmodule

module Butterfly #(parameter Jn = 1, In = 2, Bit_depth = 4) (input [Bit_depth*4-1:0] YZ, output [Bit_depth*4-1:0] XAXB);

wire [Bit_depth*2-1:0] Y;
wire [Bit_depth*2-1:0] Z;
wire [Bit_depth*2-1:0] ZW;
wire [63:0] S = 64'b0000000000_0000000000_0000000000_00__1111110000_0000000000_0000000000_00;
wire [Bit_depth*2-1:0] ZWS;
wire [Bit_depth*2-1:0] XA;
wire [Bit_depth*2-1:0] XB;
wire [63:0] W;

KN #( .Bit_depth(Bit_depth*4)) KN_1 ( .A(YZ), .B(Y), .C(Z));
Exhibitor #(.Jn(Jn), .In(In)) Exhibitor_W  ( .EXP(W));
MC #( .Bit_depth(Bit_depth)) MCZW ( .A(W), .B(Z), .C(ZW));
SUM #( .Bit_depth(Bit_depth)) SUM_XAD ( .A(ZW), .B(Y), .C(XA));
MC #( .Bit_depth(Bit_depth)) MCZWS ( .A(S), .B(ZW), .C(ZWS));
SUM #( .Bit_depth(Bit_depth)) SUM_XBD ( .A(ZWS), .B(Y), .C(XB));

assign XAXB = {XB, XA};

endmodule

module MC #(parameter Bit_depth = 4) (input [63:0] A, input [2*Bit_depth-1:0] B, output [2*Bit_depth-1:0] C);

//Complex number multiplier

reg [31:0] AR;
reg [31:0] AI;
reg [Bit_depth-1:0] BR;
reg [Bit_depth-1:0] BI;
wire [Bit_depth-1:0] CRD;
wire [Bit_depth-1:0] CID;
wire [Bit_depth-1:0] ARBR;
wire [Bit_depth-1:0] AIBI;
wire [Bit_depth-1:0] AIBIS;
wire [Bit_depth-1:0] AIBR;
wire [Bit_depth-1:0] ARBI;
wire signed [31:0] AE = 32'b1_11111_0000000000_0000000000_000000;

integer i;

always @(*) 
    begin
        for (i = 0; i<= 31; i=i+1)
        AR[i] <= A[i];   
    end

always @(*) 
    begin
        for (i = 0; i <= 31; i=i+1)
        AI[i] <= A[32+i];   
    end

always @(*) 
    begin
        for (i = 0; i<=Bit_depth-1; i=i+1)
        BR[i] <= B[i];   
    end

always @(*) 
    begin
        for (i = 0; i<=Bit_depth-1; i=i+1)
        BI[i] <= B[Bit_depth+i];   
    end



MP #( .Bit_depth(Bit_depth)) MP_ARBR ( .A(AR), .B(BR), .C(ARBR));
MP #( .Bit_depth(Bit_depth)) MP_AIBI ( .A(AI), .B(BI), .C(AIBI));
MP #( .Bit_depth(Bit_depth)) MP_AIBR ( .A(AI), .B(BR), .C(AIBR));
MP #( .Bit_depth(Bit_depth)) MP_ARBI ( .A(AR), .B(BI), .C(ARBI));
MP #( .Bit_depth(Bit_depth)) MP_AIBIS ( .A(AE), .B(AIBI), .C(AIBIS));

assign CRD = ARBR + AIBIS;
assign CID = AIBR + ARBI;

assign C = {CID,CRD};

endmodule


module MP #(parameter Bit_depth = 8) (input signed [31:0] A, input signed [Bit_depth-1:0] B, output signed [Bit_depth-1:0] C);

// Multiplier with rounding
// A(Sign+5.26)*B(Sign+Bit_depth-1)

reg [Bit_depth-1:0] Whole_part;
reg [25:0] Fractional_part;
reg [Bit_depth-1:0] Rounding_operator;
wire [Bit_depth-1:0] G = {{(Bit_depth-1){1'b0}}, {1'b1}};
wire [Bit_depth-1:0] R = {(Bit_depth){1'b0}};
wire t = Whole_part [Bit_depth-1];
wire y = Fractional_part[25];
wire [1:0] sel = {t, y};

wire [31+Bit_depth:0] Result;
integer i;

assign Result = A*B;

    always @(*) 
        begin
            for (i = 0; i <= Bit_depth-1; i=i+1)
            Whole_part[i] <= Result[26+i];   
        end

    always @(*) 
        begin
            for (i = 0; i <= 25; i=i+1)
            Fractional_part[i] <= Result[i]; 
        end


    always @(sel, G, R) begin
    case (sel)
       2'b01 : Rounding_operator = G;
       2'b11 : Rounding_operator = G;
        default: Rounding_operator = R;
    endcase 
    end

    assign C = Whole_part + Rounding_operator;

endmodule

module SUM #(parameter Bit_depth = 4) (input [Bit_depth*2-1:0] A, B, output [Bit_depth*2-1:0] C);

//Complex number adder

reg [Bit_depth-1:0] AR;
reg [Bit_depth-1:0] AI;
reg [Bit_depth-1:0] BR;
reg [Bit_depth-1:0] BI;
wire [Bit_depth-1:0] CR;
wire [Bit_depth-1:0] CI;
integer i;

always @(*) 
    begin
        for (i = 0; i<=Bit_depth-1; i=i+1)
        AR[i] <= A[i];   
    end

always @(*) 
    begin
        for (i = 0; i<=Bit_depth-1; i=i+1)
        AI[i] <= A[Bit_depth+i];   
    end

always @(*) 
    begin
        for (i = 0; i<=Bit_depth-1; i=i+1)
        BR[i] <= B[i];   
    end

always @(*) 
    begin
        for (i = 0; i<=Bit_depth-1; i=i+1)
        BI[i] <= B[Bit_depth+i];   
    end

 
assign CR = AR+BR;
assign CI = AI+BI;
assign C = {CI, CR};

endmodule

module Exhibitor #(parameter Jn = 1, In = 2) (output [63:0] EXP);

// Sign+5.26

    wire [31:0] RPE;
    wire [31:0] IPE;
    wire [31:0] Arg;
    wire [31:0] Iw;
    wire [31:0] Rw;
    wire SignSin;
    wire SignCos;
    reg [30:0] Rv;
    reg [30:0] Iv;
    wire [30:0] Ra;
    wire [30:0] Ia;
    integer i;

    exponent_argument #( .Jn(Jn), .In(In)) exponent_argument_Arg ( .Arg(Arg), .SignSin(SignSin), .SignCos(SignCos));
    Sin Sin_IPE ( .Arg(Arg), .Sin(Iw));
    Cos Cos_RPE ( .Arg(Arg), .Cos(Rw));
    
    always @(*) begin
        if (SignCos == 0) begin
                for (i = 0; i < 31; i=i+1)
                Rv[i] <= Rw[i+1];   
                end

        else begin
                for (i = 0; i < 31; i=i+1)
                Rv[i] <= ~Rw[i+1];
                end
    end
    

    assign Ra = SignCos ? Rv +1'd1 : Rv;

    always @(*) begin
        if (SignSin == 0) begin
                for (i = 0; i < 31; i=i+1)
                Iv[i] <= Iw[i+1];    
                end
                
            else begin
                for (i = 0; i < 31; i=i+1)
                Iv[i] <= ~Iw[i+1];    
                end
    end

    assign Ia = SignSin ? Iv +1'd1 : Iv;


    assign RPE = {SignCos, Ra};
    assign IPE = {SignSin, Ia};

assign EXP = {IPE, RPE};

endmodule

module exponent_argument #(parameter Jn = 1, In = 2) (output reg [31:0] Arg, output SignSin, SignCos);


              wire [31:0] Nw = (2**(Jn+1))*2;           
              wire [31:0] naw = In%(2**(Jn+1));
             
    wire [31:0] Pi = 32'b011_0010010000_1111110110_000000000;

	wire [31:0] B = Nw/4;
    wire [31:0] C = Nw/2;
	 
            function automatic [31:0] nb (input [31:0] A, R); 
                nb = ((A/R)%2)*(R-2*(A%R))+(A%R);
            endfunction

            function automatic [31:0] SC (input [31:0] J, K); 
                SC = ((J+(3*K/2))/K)%2;
            endfunction

            function automatic [31:0] SS (input [31:0] J, K); 
                SS = (J/K)%2;
            endfunction

	wire [63:0] M = 2*(Pi*nb(naw, B))/Nw;

    assign SignCos = ~SC(naw, C);
    assign SignSin = ~SS(naw, C);

	integer i;
        always @(*) 
        begin
            for (i = 0; i <= 31; i=i+1)
            Arg[i] <= M[i+2];   
        end 


endmodule

// 0000_0/00000.0000_0000000000_0000000000_000/0000000_0000000000_0000000000
//      9/87654.3210_9876543210_9876543210_987/6543210_9876543210_9876543210 


module multiplier (input [31:0] A, B, output reg [31:0] C);

// 5.27

wire [63:0] Result;
integer i;

assign Result = A*B;

    always @(*) 
    begin
        for (i = 0; i<=31; i=i+1)
        C[i] <= Result[27+i];   
    end

endmodule

module Sin (input [31:0] Arg, output [31:0] Sin);
// 1/3! (27 characters after the dot)
wire [31:0] CA = 32'b00000_0010101010_1010101010_1010101;
// 1/5! (27 characters after the dot)
wire [31:0] CB = 32'b00000_0000001000_1000100010_0010001;
// 1/7! (27 characters after the dot)
wire [31:0] CC = 32'b00000_0000000000_0011010000_0000110;
wire [31:0]AM1;
wire [31:0]AM2;
wire [31:0]AM3;
wire [31:0]AM4;
wire [31:0]AM5;
wire [31:0]AM6;
wire [31:0]AM7;
assign AM1 = Arg;
multiplier multiplier_AM2 ( .A(AM1), .B(AM1), .C(AM2));
multiplier multiplier_AM3 ( .A(AM1), .B(AM2), .C(AM3));
multiplier multiplier_AM4 ( .A(AM1), .B(AM3), .C(AM4));
multiplier multiplier_AM5 ( .A(AM1), .B(AM4), .C(AM5));
multiplier multiplier_AM6 ( .A(AM1), .B(AM5), .C(AM6));
multiplier multiplier_AM7 ( .A(AM1), .B(AM6), .C(AM7));
wire [31:0]A;
wire [31:0]B;
wire [31:0]C;
multiplier multiplier_A ( .A(AM3), .B(CA), .C(A));
multiplier multiplier_B ( .A(AM5), .B(CB), .C(B));
multiplier multiplier_C ( .A(AM7), .B(CC), .C(C));

assign Sin = (AM1-A)+(B-C);

endmodule

module Cos (input [31:0] Arg, output [31:0] Cos);
// 1/2! (27 characters after the dot)
wire [31:0] CA = 32'b00000_1000000000_0000000000_0000000;
// 1/4! (27 characters after the dot)
wire [31:0] CB = 32'b00000_0000101010_1010101010_0101010;
// 1/6! (27 characters after the dot)
wire [31:0] CC = 32'b00000_0000000001_0110110000_0101101;
// 1/8! (27 characters after the dot)
//wire [31:0] CD = 32'b00000_0000000000_0000011010_0000000;
wire [31:0] CD = 32'b00000_0000000000_101110011_111111;
wire [31:0] Const = 32'b00001_0000000000_0000000000_0000000;
wire [31:0]AM1 = Arg;
wire [31:0]AM2;
wire [31:0]AM3;
wire [31:0]AM4;
wire [31:0]AM5;
wire [31:0]AM6;
wire [31:0]AM7;
wire [31:0]AM8;
assign AM1 = Arg;
multiplier multiplier_AM2 ( .A(AM1), .B(AM1), .C(AM2));
multiplier multiplier_AM3 ( .A(AM1), .B(AM2), .C(AM3));
multiplier multiplier_AM4 ( .A(AM1), .B(AM3), .C(AM4));
multiplier multiplier_AM5 ( .A(AM1), .B(AM4), .C(AM5));
multiplier multiplier_AM6 ( .A(AM1), .B(AM5), .C(AM6));
multiplier multiplier_AM7 ( .A(AM1), .B(AM6), .C(AM7));
multiplier multiplier_AM8 ( .A(AM1), .B(AM7), .C(AM8));
wire [31:0]A;
wire [31:0]B;
wire [31:0]C;
wire [31:0]D;
multiplier multiplier_A ( .A(AM2), .B(CA), .C(A));
multiplier multiplier_B ( .A(AM4), .B(CB), .C(B));
multiplier multiplier_C ( .A(AM6), .B(CC), .C(C));
multiplier multiplier_D ( .A(AM8), .B(CD), .C(D));

assign Cos = Const+D+B-A-C;

endmodule

module KN #(parameter Bit_depth = 8) (input [Bit_depth-1:0]A, output reg [(Bit_depth/2)-1:0] B, C);

integer i;

always @(*) 
    begin          
        for (i = 0; i<=(Bit_depth/2)-1; i=i+1)
        B[i] <= A[i];                   
    end

always @(*) 
    begin          
        for (i = 0; i<=(Bit_depth/2)-1; i=i+1)
        C[i] <= A[i+(Bit_depth/2)];                  
    end

endmodule

module  R #(parameter N) (input [N-1:0] A, output [N-1:0] B);

assign B = A;

endmodule

module Permutation_conductor #(parameter Bus = 8, Bit_depth = 6, Jn = 0) (input [Bus*Bit_depth-1:0] A, output reg [Bus*Bit_depth-1:0] B);

function [$clog2(Bus)-1:0] permutation (input [$clog2(Bus)-1:0] I);

integer Jt;
integer k;
reg [$clog2(Bus)-1:0] A;
reg [$clog2(Bus)-1:0] B; 
reg [$clog2(Bus)-1:0] C; 
reg [$clog2(Bus)-1:0] D; 
reg [$clog2(Bus)-1:0] Q; 
reg [$clog2(Bus)-1:0] W; 
reg [$clog2(Bus)-1:0] E; 
reg [$clog2(Bus)-1:0] R; 
reg [$clog2(Bus)-1:0] T; 

begin

Jt = 2**(Jn+2);
R = 0;

for ( k = 0 ; k <= I ; k=k+1 ) begin

A = k%Jt;
if(A == 0) B = Jt/2; else B = 0;
if(k == 0) C = 0; else C = 1;
D =B*C;
Q = k%2;
W = (k+1)%2;
E = C*W;
R = E+D+R;
T = R + Jt/2;
end

permutation = W*R+T*Q;

end

endfunction


integer i;
integer j;

always @(*) begin
    for (i = 0 ; i < Bus; i = i + 1) begin
        for (j = 0; j < Bit_depth ; j = j + 1) begin
           B[i*Bit_depth+j] = A[permutation(i)*Bit_depth+j];
        end
    end  
end

endmodule

module Permutation_conductor_two #(parameter Bus = 8, Bit_depth = 4, Jn = 1) (input [Bus*Bit_depth-1:0] A, output reg [Bus*Bit_depth-1:0] B);

function [$clog2(Bus)-1:0] permutation (input [$clog2(Bus)-1:0] I);

integer Jt; 

begin

Jt = 2**Jn;
if ((((I+Jt)/Jt)+2)%2 == 0) permutation = I - (Jt-1-(I+Jt)%Jt); 
else permutation = I+(I+Jt)%Jt;

end

endfunction

integer i;
integer j;

        always @(*) 
        begin
            for (i = 0 ; i < Bus; i = i + 1) 
            begin
                for (j = 0; j < Bit_depth ; j = j + 1) 
                    begin
                        B[i*Bit_depth+j] = A[permutation(i)*Bit_depth+j];
                    end
            end 
        end

endmodule