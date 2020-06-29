classdef MarchingCube
    properties
orient;
cell;
    end
    
    methods
        function obj = MarchingCube( )
            
            for i=0:255
                str = dec2bin(i,8);
                A =str(0);
                B =str(1);
                C =str(2);
                D =str(3);

                E =str(4);
                F =str(5);
                G =str(6);
                H =str(7);
                
                baseType[i] = getMatch(A,B,C,D,E,F,G,H); 
                if(baseType!=-1) 
                    orient[i] = [0;0;0];
                end
                
                    %X=0,Y=0,Z=0.
                    if(Xorient(A,B,C,D,E,F,G,H,i,1,0,0)) continue;end; 
                    if(Xorient(A,B,C,D,E,F,G,H,i,2,0,0)) continue;end;
                    if(Xorient(A,B,C,D,E,F,G,H,i,3,0,0)) continue;end;
                    if(Xorient(A,B,C,D,E,F,G,H,i,0,0,0)) continue;end;

                %X=0,Y=0,Z=0.
                if(Yorient(A,B,C,D,E,F,G,H,i,0,1,0)) continue;end; 
                %X=0,Y=1,Z=0.

                    if(Zorient(A,B,C,D,E,F,G,H,i,0,1,1)) continue;end; 
                    if(Zorient(A,B,C,D,E,F,G,H,i,0,1,2)) continue;end;
                    if(Zorient(A,B,C,D,E,F,G,H,i,0,1,3)) continue;end;
                    if(Zorient(A,B,C,D,E,F,G,H,i,0,1,4)) continue;end;

                %X=0,Y=1,Z=0.
                if(Yorient(A,B,C,D,E,F,G,H,i,0,2,0)) continue;end; 
                %X=0,Y=2,Z=0.

                    if(Xorient(A,B,C,D,E,F,G,H,i,1,2,0)) continue;end; 
                    if(Xorient(A,B,C,D,E,F,G,H,i,2,2,0)) continue;end;
                    if(Xorient(A,B,C,D,E,F,G,H,i,3,2,0)) continue;end;
                    if(Xorient(A,B,C,D,E,F,G,H,i,0,2,0)) continue;end;

                %X=0,Y=2,Z=0.
                if(Yorient(A,B,C,D,E,F,G,H,i,0,3,0)) continue;end; 
                %X=0,Y=3,Z=0.
                
                    if(Zorient(A,B,C,D,E,F,G,H,i,0,3,1)) continue;end; 
                    if(Zorient(A,B,C,D,E,F,G,H,i,0,3,2)) continue;end;
                    if(Zorient(A,B,C,D,E,F,G,H,i,0,3,3)) continue;end;
                    if(Zorient(A,B,C,D,E,F,G,H,i,0,3,0)) continue;end;
                
                %X=0,Y=3,Z=0.
                if(Yorient(A,B,C,D,E,F,G,H,i,0,0,0)) continue;end; 
                %X=0,Y=0,Z=0.

                %X=0,Y=0,Z=0.
                if(Zorient(A,B,C,D,E,F,G,H,i,0,0,1)) continue;end;
                %X=0,Y=0,Z=1.

                    if(Yorient(A,B,C,D,E,F,G,H,i,0,1,1)) continue;end; 
                    if(Yorient(A,B,C,D,E,F,G,H,i,0,2,1)) continue;end; 
                    if(Yorient(A,B,C,D,E,F,G,H,i,0,3,1)) continue;end; 
                    if(Yorient(A,B,C,D,E,F,G,H,i,0,0,1)) continue;end; 


                %X=0,Y=0,Z=1.
                if(Zorient(A,B,C,D,E,F,G,H,i,0,0,2)) continue;end;
                if(Zorient(A,B,C,D,E,F,G,H,i,0,0,3)) continue;end;
                %X=0,Y=0,Z=3.
                
                    if(Yorient(A,B,C,D,E,F,G,H,i,0,1,3)) continue;end; 
                    if(Yorient(A,B,C,D,E,F,G,H,i,0,2,3)) continue;end; 
                    if(Yorient(A,B,C,D,E,F,G,H,i,0,3,3)) continue;end; 
                    if(Yorient(A,B,C,D,E,F,G,H,i,0,0,3)) continue;end; 
    
                %X=0,Y=0,Z=3.
                if(Zorient(A,B,C,D,E,F,G,H,i,0,0,0)) continue;end;
                %X=0,Y=0,Z=0.
                
            end
        end
        function getMatch(A,B,C,D,E,F,G,H) 
            cellType = A*128 + B * 64 + C*32+ D*16 + E*8 + F*4 + G*2 + F;   

            if(cellType==13) return 5; end
            if(cellType==13) return 5; end
            if(cellType==13) return 5; end
            if(cellType==13) return 5; end
        end
        
        function fourSwap(i,j,k,l)
            temp=l;
            l=k;
            k=j;
            j=i;
            i=temp;
        end
        
        function found = Xorient(A,B,C,D,E,F,G,H,i,x,y,z)
            fourSwap(A,B,C,D);
            fourSwap(E,F,G,H);

            baseType[i] = getMatch(A,B,C,D,E,F,G,H); 
            if(baseType!=-1) 
                orient[i] = [x;y;z];
                return 1;
            end
            return 0;
        
        end
        function Yorient(A,B,C,D,E,F,G,H,i,x,y,z)
            fourSwap(A,B,C,D);
            fourSwap(E,F,G,H);

            baseType[i] = getMatch(A,B,C,D,E,F,G,H); 
            if(baseType!=-1) 
                orient[i] = [x;y;z];
                return 1;
            end
            return 0;

            
        end
        
        function Zorient(A,B,C,D,E,F,G,H,i,x,y,z)
            fourSwap(A,B,C,D);
            fourSwap(E,F,G,H);

            baseType[i] = getMatch(A,B,C,D,E,F,G,H); 
            if(baseType!=-1) 
                orient[i] = [x;y;z];
                return 1;
            end
            return 0;
            
        end
        
        function getBasetype( celltype)
            
            return baseType(celltype)
        end
        
        function getOrientation(celltype)
            return orient(cellType)
        end
    
    end
end

