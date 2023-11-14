clear all; close all;

columNum=8; % 8

fileNameHat=['relax'];  % relax Gr40x3_4L2_CNT5@5x2_PBD_NPTRlx
% fileNameHat=['holdAF_'];  % _NPTCG relax100000
fileNameGroup=[65000]; %   Gr40x6_7L2_CNT13@13x2_PBD-MinCG

nFile=size(fileNameGroup,2);
for  i=1:nFile
    iFileName=[fileNameHat num2str(fileNameGroup(1,i))];
    fre01 = fopen([ iFileName '.dump'],'r');
    tline = fgetl(fre01);
    tline = fgetl(fre01);
    tline = fgetl(fre01);
    tline = fgetl(fre01);
    nAtom=str2double(tline);
    
    tline = fgetl(fre01);  % ----------
    
    tline = fgetl(fre01);
    iSizeX=sscanf(tline,'%f')';
    tline = fgetl(fre01);
    iSizeY=sscanf(tline,'%f')';
    tline = fgetl(fre01);
    iSizeZ=sscanf(tline,'%f')';
    
    boxsize=[iSizeX; iSizeY; iSizeZ];
    
    iCoorTemp=zeros(nAtom,columNum);
    
    tline = fgetl(fre01);  % ----------
    for j=1:nAtom
        tline = fgetl(fre01);
        iLine=sscanf(tline,'%f')';
        
        iCoorTemp(j,:)=iLine;

    end
    iCoor=sortrows(iCoorTemp,1);
    
    atomType=max(iCoor(:,3));
    
    writeData(nAtom,iCoor,boxsize,iFileName,atomType)
    
    disp('Done');

end

fclose(fre01); 



function []=writeData(nAtom,iCoor,boxsize,iFileName,atomType)

nAtom_full=nAtom;
coor_full=iCoor;

% *************************************************************************
dX=0.0;
xlo=boxsize(1,1)-dX; xhi= boxsize(1,2)+dX;
ylo=boxsize(2,1)-dX; yhi= boxsize(2,2)+dX;
zlo=boxsize(3,1)-dX; zhi= boxsize(3,2)+dX;

% xlo=min(coor(:,1))-dX; xhi= max(coor(:,1))+dX;
% ylo=min(coor(:,2))-dX; yhi= max(coor(:,2))+dX;
% zlo=min(coor(:,3))-dX; zhi= max(coor(:,3))+dX;

% 
% dataFileName='Onion_L3_60C1_240C1_540C1';
fileID = fopen([iFileName '.data'],'w');
fprintf(fileID,'LAMMPS Data File. atom_stye full by Mingda\n');
fprintf(fileID,'%d atoms\n',nAtom_full);
fprintf(fileID,'0 bonds\n');
fprintf(fileID,'0 angles\n');
fprintf(fileID,'0 dihedrals\n');
fprintf(fileID,'0 impropers\n');
fprintf(fileID,'%d atom types\n',atomType);
fprintf(fileID,'0 bond types\n');
fprintf(fileID,'0 angle types\n');
fprintf(fileID,'0 dihedral types\n');
fprintf(fileID,'0 improper types\n');
fprintf(fileID,'%.4f %.4f xlo xhi\n',xlo,xhi);
fprintf(fileID,'%.4f %.4f ylo yhi\n',ylo,yhi);
fprintf(fileID,'%.4f %.4f zlo zhi\n',zlo,zhi);
fprintf(fileID,'\n');
fprintf(fileID,' Masses \n');
fprintf(fileID,'\n');

for i=1:atomType
    fprintf(fileID,'%d 12.010700 # C\n',i);
end

fprintf(fileID,'\n');
fprintf(fileID,' Atoms # full\n');
fprintf(fileID,'\n');
for i=1:nAtom_full
    fprintf(fileID,'%d %d %d 0.0000 %.4f %.4f %.4f\n',...
        i,coor_full(i,2),iCoor(i,3),iCoor(i,4),iCoor(i,5),coor_full(i,6));
end
fclose(fileID);

end