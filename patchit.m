function patchit(I,PatchSize,varargin)

%     patchIT Tool
%     A Multipurpose Patch Creation Tool for  Image Processing Applications
%     v1.1
%     
%     Copyright (C) 2022 Anastasios Kesidis akesidis@uniwa.gr
%      
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.


% --- load packages (only in Octave)
if exist('OCTAVE_VERSION', 'builtin')>0
    pkg load image;
end

dbstop if error; % <--- TEMP
dbstop if warning; % <--- TEMP

s=fprintf('Initializing...\n');


p=inputParser;
addRequired(p,'I',@(x) ischar(x) || (all(isfinite(x(:))) && all(isnumeric(x(:))) && ~isvector(x)));
addRequired(p,'PatchSize',@(x) all(isfinite(x)) && all(x==round(x)) && numel(x)==2 && all(x>0) && isrow(x));

DefaultMode='sliding';
addParameter(p,'Mode',DefaultMode,@(x) any(strcmpi(x,{'sliding','random'})));
DefaultMask=[];
addParameter(p,'Mask',DefaultMask);
DefaultCount=10;
addParameter(p,'Count',DefaultCount,@(x) isfinite(x) && x==round(x) && isscalar(x) && x>0);
DefaultStride=[0 0];
addParameter(p,'Stride',DefaultStride,@(x) all(isfinite(x)) && all(x==round(x)) && numel(x)==2 && all(x>0) && isrow(x));
DefaultPatchIndexing='default';
addParameter(p,'PatchIndexing',DefaultPatchIndexing,@(x) ischar(x));

DefaultRandomTrials=1;
addParameter(p,'RandomMaxOverlap',DefaultRandomTrials,@(x) isfinite(x) && isscalar(x) && x>=0);
DefaultRandomAttempts=100;
addParameter(p,'RandomAttempts',DefaultRandomAttempts,@(x) isfinite(x) && isscalar(x) && x>0);

DefaultOrder='row';
addParameter(p,'Order',DefaultOrder,@(x) any(strcmpi(x,{'Column','Row','ColumnRev','RowRev'})));


DefaultSaveMode='mat';
addParameter(p,'SaveMode',DefaultSaveMode,@(x) any(strcmpi(x,{'mat','images','raw'})));

DefaultSaveMatFilename='patches';
addParameter(p,'SaveMatFilename',DefaultSaveMatFilename,@(x) ischar(x));

DefaultSaveImagesTemplate='patch0';
addParameter(p,'SaveImagesTemplate',DefaultSaveImagesTemplate,@(x) ischar(x));
DefaultSaveImagesExt='png';
addParameter(p,'SaveImagesExt',DefaultSaveImagesExt,@(x) ischar(x));

DefaultSaveRawFilename='raw';
addParameter(p,'SaveRawFilename',DefaultSaveRawFilename,@(x) ischar(x));

DefaultSavePatchesPosFilename='pos';
addParameter(p,'SavePatchesPosFilename',DefaultSavePatchesPosFilename,@(x) ischar(x));

DefaultProcessingMode='direct';
addParameter(p,'ProcessingMode',DefaultProcessingMode,@(x) any(strcmpi(x,{'direct','block'})));


DefaultShow=false;
addParameter(p,'Show',DefaultShow,@(x) islogical(x));
DefaultTileRows=5;
addParameter(p,'TileRows',DefaultTileRows,@(x) isnumeric(x) && isscalar(x) && (x>0));
DefaultTileColumns=9;
addParameter(p,'TileColumns',DefaultTileColumns,@(x) isnumeric(x) && isscalar(x) && (x>0));
DefaultShowPatchesID=[];
addParameter(p,'ShowPatchesID',DefaultShowPatchesID,@(x) isnumeric(x) && isscalar(x) && (x>0));

try
    parse(p,I,PatchSize,varargin{:});

    %% Preparation
    % --- Get image dimensions and number of bands
    if ischar(I)
        if ~isfile(I)
            error('File %s not found\n',I);
        else
            I=imread(I);
        end
    end

    if ischar(p.Results.Mask)
        if ~isfile(p.Results.Mask)
            error('File %s not found\n',p.Results.Mask);
        else
            M=imread(p.Results.Mask);
            if size(M,3)>1
                M=rgb2gray(M);
            end
        end
    elseif isempty(p.Results.Mask)
        M=(ones(size(I,1),size(I,2)))==1;
    else
        M=p.Results.Mask;
    end
    if ~isequal(size(M,1:2),size(I,1:2))
        error('Invalid input parameter.\nImage and mask dimensions are incompatible.\n');
    end

    % --- Get image dimensions and number of channels
    [Ir,Ic,Ib]=size(I);

    % --- Get patch dimensions
    Pr=p.Results.PatchSize(2);
    Pc=p.Results.PatchSize(1);

    % --- Ensure patch size is not greater than image size
    if or(Ir<Pr,Ic<Pc)
        error('Invalid input parameter.\nPatch dimensions should not exceed image dimensions.\n');
    end

    % --- Define default stride equal to patch size
    if ~any(p.Results.Stride==0)
        Stride=p.Results.Stride;
    else
        Stride=PatchSize;
    end

    %% Define patches
    % --- Define the candidate patch positions determined by the Mask
    PatchMaskPool=imerode(double(padarray(M,[1 1])),ones(Pr,Pc));
    PatchMaskPool=PatchMaskPool(2:end-1,2:end-1);
    % --- Mask candidate positions
    Idx2=find(PatchMaskPool);

    if strcmp(p.Results.Mode,'sliding')
        % --- Define the starting and ending positions for the sliding parameters
        Startx=ceil(Pc/2);
        Endx=Ic-floor(Pc/2);
        Starty=ceil(Pr/2);
        Endy=Ir-floor(Pr/2);
        % --- Define SLIDING patches
        x1=Startx:Stride(1):Endx;
        y1=Starty:Stride(2):Endy;
        % --- Sliding candidate positions
        [X1,Y1]=meshgrid(x1,y1);
        Idx1=sub2ind(size(M),Y1(:),X1(:));
        % --- Final patch center positions
        PatchCenters=intersect(Idx1,Idx2);
    elseif strcmpi(p.Results.Mode,'random')
        % --- Define RANDOM patch locations, i.e. (x,y) coordinates of their upper left corner
        % x1=floor(rand(p.Results.PatchCount)*(C-p.Results.PatchSize(1)))+1;
        % y1=floor(rand(p.Results.PatchCount)*(R-p.Results.PatchSize(2)))+1;
        % --- Define the starting and ending positions for the sliding parameters
        Startx=ceil(Pc/2);
        Endx=Ic-floor(Pc/2);
        Starty=ceil(Pr/2);
        Endy=Ir-floor(Pr/2);
        % --- Prepare RANDOM candidate patches
        x1=Startx:Endx;
        y1=Starty:Endy;
        [X1,Y1]=meshgrid(x1,y1);
        Idx1=sub2ind(size(M),Y1(:),X1(:));
        Weights=double(M(Idx1));
        % --- Find the candidate patches
        for i=1:p.Results.RandomAttempts
            PatchCenters=datasample(Idx1,p.Results.Count,'weights',Weights);

            % --- Prepare the 4 corners for all patches
            [r,c]=ind2sub([Ir Ic],PatchCenters);
            x1=c-ceil(Pc/2)+1;
            y1=r-ceil(Pr/2)+1;
            x2=x1+Pc-1;
            y2=y1+Pr-1;
            A=[x1 y1 (x2-x1+1) (y2-y1+1)];
            Overlap=rectint(A,A);
            Overlap=Overlap.*~eye(size(Overlap));
            OverlapFound=any(Overlap(:)>prod(p.Results.PatchSize)*p.Results.RandomMaxOverlap);
            if ~OverlapFound
                break;
            end
        end
        if OverlapFound
            error('The maximum number of attempts for random patches is reached .\n');
        end
    end

    % --- Final number of patches
    N=length(PatchCenters);

    %% Prepare indexing
    UpLeftCorners=uint32(PatchCenters-(ceil(Pr/2)-1)-Ir*(ceil(Pc/2)-1));
    [r,c]=ind2sub([Pr Pc],patch_indexing(Pr,Pc,p.Results.PatchIndexing));
    PatchIdx=uint32(sub2ind([Ir Ic],r,c));
    UpLeftCornersStep=length(UpLeftCorners)+100;
    if strcmpi(p.Results.ProcessingMode,'block')
        UpLeftCornersStep=ceil(sqrt(length(UpLeftCorners)));
    end
    BandsOffset=uint32((0:Ib-1)*Ir*Ic);

    %% Prepare order
    [r,c]=ind2sub([Ir Ic],UpLeftCorners);
    if strcmpi(p.Results.Order,'Row')
        tmp=r;
        r=c;
        c=tmp;
    elseif strcmpi(p.Results.Order,'ColumnRev')
        r=flipud(r);
        c=flipud(c);
    elseif strcmpi(p.Results.Order,'RowRev')
        tmp=r;
        r=c;
        c=tmp;
        r=flipud(r);
        c=flipud(c);        
    end
    UpLeftCorners=uint32(sub2ind([Ir Ic],r,c));

    %% Prepare save images
    if strcmpi(p.Results.SaveMode,'images')
        % --- Define save folder and filenames
        [SaveImagesPath,~,SaveImagesExt]=fileparts(p.Results.SaveImagesTemplate);
        if isempty(SaveImagesExt)
            SaveImagesExt=['.' p.Results.SaveImagesExt];
        end
        [reg_s1,reg_s2]=regexp(p.Results.SaveImagesTemplate, '0{1,}','match','split');
        if ~isempty(reg_s1)
            PatchSavePattern=['%s%0' num2str(length(reg_s1{:})) 'd%s'];
        else
            PatchSavePattern='%s%1d%s';
        end
        % --- Ensure that destination folder exists
        [~,~]=mkdir(SaveImagesPath);
        % --- Delete any existing patch files
        delete(sprintf('%s%s%s',reg_s2{1},'*',SaveImagesExt));
    end

    %% Prepare save mat
    if strcmpi(p.Results.SaveMode,'mat')
        [SaveMatPath,SaveMatFilename,SaveMatExt]=fileparts(p.Results.SaveMatFilename);
        if isempty(SaveMatExt)
            SaveMatExt='.mat';
        end
        if isempty(SaveMatPath)
            MatFile=[SaveMatFilename,SaveMatExt];
        else
            MatFile=[SaveMatPath,filesep,SaveMatFilename,SaveMatExt];
        end
        % --- Ensure that destination folder exists
        [~,~]=mkdir(SaveMatPath);
        % --- Delete any existing mat file
        if exist(MatFile,'file')
            delete(MatFile);
        end
    end
    %% Prepare save raw
    if strcmpi(p.Results.SaveMode,'raw')
        [SaveRawPath,SaveRawFilename,SaveRawExt]=fileparts(p.Results.SaveRawFilename);
        if isempty(SaveRawExt)
            SaveRawExt='.txt';
        end
        if isempty(SaveRawPath)
            FileSeparator='';
        else
            FileSeparator=filesep;
        end


        if Ib==1
            RawFile{1}=[SaveRawPath,FileSeparator,SaveRawFilename,SaveRawExt];
            if exist(RawFile{1},'file')
                delete(RawFile{1});
            end
        else
            % --- Ensure that destination folder exists
            [~,~]=mkdir(SaveRawPath);
            for i=1:Ib
                RawFile{i}=[SaveRawPath,FileSeparator,[SaveRawFilename sprintf(['%0' num2str(length(num2str(Ib))) 'd'],i)],SaveRawExt];
                if exist(RawFile{i},'file')
                    delete(RawFile{i});
                end
            end
        end
    end

    %% Prepare save patches pos
    % --- in all save cases a pos.txt file is created
    [SavePatchesPosPath,SavePatchesPosFilename,SavePatchesPosExt]=fileparts(p.Results.SavePatchesPosFilename);
    if isempty(SavePatchesPosExt)
        SavePatchesPosExt='.txt';
    end
    PatchesPosFile=[SavePatchesPosPath,SavePatchesPosFilename,SavePatchesPosExt];
    if exist(PatchesPosFile,'file')
        delete(PatchesPosFile);
    end

    fprintf(repmat('\b',1,s));
    fprintf('Initializing OK\n');

    %% Save main

    % --- Initialize
    k=0;
    s=0;
    c=0;
    if UpLeftCornersStep>length(UpLeftCorners)
        s=fprintf('Saving %g patches...',N);
    end
    tic;

    for i=1:UpLeftCornersStep:length(UpLeftCorners)
        c=c+1;
        StartID=i;
        EndID=min(i+UpLeftCornersStep-1,length(UpLeftCorners));
        UpLeftCornersPart=UpLeftCorners(StartID:EndID);
        Idx1=bsxfun(@plus,PatchIdx,(UpLeftCornersPart-1)');
        Idx3=bsxfun(@plus,Idx1,reshape(BandsOffset,1,1,[]));
        np=length(UpLeftCornersPart);
        ImagePatchesIdx=reshape(pagetranspose(permute(Idx3,[3 1 2])),Pr,Pc,[],np);

        [y1,x1]=ind2sub([Ir,Ic],UpLeftCornersPart);
        x2=x1+Pc-1;
        y2=y1+Pr-1;
        PatchPos=[x1 y1 x2 y2];

        if strcmpi(p.Results.SaveMode,'images') % --- Save images
            for j=StartID:EndID
                k=k+1;
                PatchSaveFilename=sprintf(PatchSavePattern,reg_s2{1},k,SaveImagesExt);
                imwrite(I(ImagePatchesIdx(:,:,:,j-StartID+1)),PatchSaveFilename);
            end

        elseif strcmpi(p.Results.SaveMode,'mat') % --- Save mat
            if ~exist('m','var')
                patch(:,:,:,StartID:EndID)=I(ImagePatchesIdx);
                pos(StartID:EndID,:)=PatchPos;
                save(MatFile,'patch','-v7.3');
                m = matfile(MatFile,'Writable',true);
            else
                m.patch(:,:,:,StartID:EndID)=I(ImagePatchesIdx);
                % --- NOTICE !! you can read PART of a HUGE mat file as follows:
                % --- m = matfile('patches.mat');
                % --- y=m.patch(:,:,:,1:100); <--- reads only the first 100 patches
            end

        elseif strcmpi(p.Results.SaveMode,'raw') % --- Save raw
            % --- Separate text file for each band
            % --- Each patch is a line of the text file
            for j=1:Ib
                fid=fopen(RawFile{j},'at');
                fprintf(fid,[repmat('%3d ',1,size(Idx3,1)) '\n'],I(Idx3(:,:,j)));
                fclose(fid);
            end

        end
        % --- Save patches position
        % --- Each line of the text file is the patch coordinates [x1 y1 x2 y2]
        fid=fopen(PatchesPosFile,'at');
        fprintf(fid,[repmat('%3d ',1,4) '\n'],PatchPos');
        fclose(fid);

        % --- update message in Command Window
        if mod(c,10)==0
            fprintf(repmat('\b',1,s));
            t=toc;
            s=fprintf('Saving %g/%g patches. Estimated remaining time: %0.1f seconds...\n',i+np-1,N,t/(i+np-1)*(N-(i+np-1)));
        end
    end
    fprintf(repmat('\b',1,s));
    fprintf('Saving %g patches OK\n',N);

catch ME
    fprintf('\n*** PATCHIT message ***\n');
    switch ME.identifier
        case {'MATLAB:nomem','MATLAB:array:SizeLimitExceeded'}
            fprintf('The requested memory is too large.\nConsider setting parameter ''ProcessingMode'' to ''Block''.\n');
        case {'MATLAB:save:noParentDir'}
            fprintf('The output mat file\n%s\ncannot be created.\n',MatFile);
        otherwise
            if ~isempty(ME.message)
                fprintf(ME.message);
            else
                fprintf('Unknown error\n');
            end
    end
    return;
end



function f=patch_indexing(Pr,Pc,Mode)

P=reshape(1:Pr*Pc,Pr,Pc);

ModeParts=regexp(Mode,'\w+','match');
flag=0;

while ~isempty(ModeParts)
    Part=ModeParts{1};
    ModeParts=ModeParts(2:end);

    if strcmpi(Part,'default')
    elseif strcmpi(Part,'swap')
        P=P';
    elseif strcmpi(Part,'mirror')
        P=fliplr(P);
    elseif strcmpi(Part,'flip')
        P=flipud(P);
    elseif strcmpi(Part,'spiralout')
        if Pr~=Pc
            error('Spiral patch indexing can be applied only on square patches.\n');
        end
        P=spiral(Pc);
        flag=1;
    elseif strcmpi(Part,'spiralin')
        if Pr~=Pc
            error('Spiral patch indexing can be applied only on square patches.\n');
        end
        P=fliplr(Pr*Pc+1-spiral(Pc));
        flag=1;
    else
        error('Invalid patch indexing mode.\n');
    end
end
if flag==0
    f=P(:);
else
    [~,f]=sort(P(:));
end


