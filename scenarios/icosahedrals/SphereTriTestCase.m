classdef SphereTriTestCase <matlab.unittest.TestCase
    % $Author: Peter Gagarinov, PhD  <pgagarinov@gmail.com> $
    % $Copyright: Peter Gagarinov, PhD,
    %            Moscow State University,
    %            Faculty of Computational Mathematics and Computer Science,
    %            System Analysis Department 2011-2016 $
    %
    properties(Access=private)
        rootDataDir
    end
    properties (Constant, GetAccess=private)
        TRI1_VERT=[1 0 0;0 1 0;0 0 1];
        TRI1_FACE=[1 2 3];
        %
        TRI2_VERT=[1 0 0;0 1 0;0 0 1;1 0 1];
        TRI2_FACE=[1 2 3;1 3 4;1 2 4;2 3 4];
        %
        TRI3_VERT=[0 0 0;1 0 0;0 1 0;-0.5 0 0;0 -0.5 0];
        TRI3_FACE=[1 2 3;1 3 4;1 5 2];
        TRI3_EDGE=[1 2;1 3;1 4;1 5;2 3;2 5;3 4];
        TRI3_F2E=[1 5 2;2 7 3;4 6 1];
        TRI3_F2E_DIR=[true true true;true true true;true false true];
        
        %
        TRI31_VERT=[0 0 0;1 0 0;0 1 0;-0.5 0 0];
        TRI31_FACE=[1 2 3;1 3 4];
    end
    methods
        function dirName=getDataDir(self)
            dirName=self.rootDataDir;
        end
        function S=loadData(self,fileName)
            S=load([self.getDataDir,filesep,fileName]);
        end
        %
        function saveData(self,fileName,SInp) %#ok<INUSD>
            save([self.getDataDir,filesep,fileName],'-struct','SInp');
        end
        %
        function self = SphereTriTestCase(varargin)
            self = self@matlab.unittest.TestCase(varargin{:});
            self.rootDataDir=[fileparts(mfilename('fullpath')),...
                filesep,'TestData'];
        end
        function [vMat,fMat,SStats]=aux_shrinkfacetri(self,vMat,fMat,varargin)
            import gras.geom.tri.*;
            [vMat,fMat,SStats,eMat,f2eMat]=shrinkfacetri(...
                vMat,fMat,varargin{:});
            nEdges=size(eMat,1);
            self.assertEqual(nEdges,length(unique(f2eMat)));
        end
        %
    end
    methods (Test)
        function testShrinkFaceTriOneFace(self)
            vMat=self.TRI1_VERT;
            fMat=self.TRI1_FACE;
            [v1Mat,f1Mat]=self.aux_shrinkfacetri(vMat,fMat,0,1);
            v1ExpMat=[1 0 0;0 1 0;0 0 1;0.5 0.5 0;0.5 0 0.5;0 0.5 0.5];
            f1ExpMat=[4 6 5;4 5 1;6 4 2;5 6 3];
            self.assertTrue(isequal(v1Mat,v1ExpMat));
            self.assertTrue(isequal(f1Mat,f1ExpMat));
            fMat=[3 2 1];
            [~,~]=self.aux_shrinkfacetri(vMat,fMat,0,1);
            fMat=[2 3 1];
            [~,~]=self.aux_shrinkfacetri(vMat,fMat,0,1);
            fMat=[3 1 2];
            [~,~]=self.aux_shrinkfacetri(vMat,fMat,0,1);
            [~,~]=self.aux_shrinkfacetri(vMat,fMat,0,4);
            %
        end
        function testShrinkFaceTri3Faces(self)
            vMat=self.TRI2_VERT;
            fMat=self.TRI2_FACE;
            [~,~]=self.aux_shrinkfacetri(vMat,fMat,0,2);
        end
        %
        function testShrinkFaceTri2Face1Part(self)
            vMat=self.TRI31_VERT;
            fMat=self.TRI31_FACE;
            [~,~]=self.aux_shrinkfacetri(vMat,fMat,realsqrt(2)-0.001,1);
        end
        %
        function testShrinkFaceTri3Face1Part(self)
            import gras.geom.tri.*;
            vMat=self.TRI3_VERT;
            fMat=self.TRI3_FACE;
            [v1Mat,f1Mat]=self.aux_shrinkfacetri(vMat,fMat,realsqrt(2)-0.001,...
                1,@(x)(x+repmat([0 0 0.2],size(x,1),1)));
            %
            isFaceThereVec=isface(v1Mat,f1Mat,[1 6 2;1 7 3]);
            self.assertTrue(all(isFaceThereVec));
            %
            [~,~]=self.aux_shrinkfacetri(v1Mat,f1Mat,0,...
                3,@(x)(x+repmat([0 0 0.2],size(x,1),1)));
        end
        %
        function testIsFace(self)
            import gras.geom.tri.*;
            %
            vMat=self.TRI3_VERT;
            fMat=self.TRI3_FACE;
            isFaceThereVec=isface(vMat,fMat,[1 5 4;fMat;2 5 4]);
            self.assertTrue(isequal(isFaceThereVec,...
                [false;true;true;true;false]));
        end
        %
        function testMapFace2Edge(self)
            import gras.geom.tri.*;
            vMat=self.TRI3_VERT;
            fMat=self.TRI3_FACE;
            eExpMat=self.TRI3_EDGE;
            f2eExpMat=self.TRI3_F2E;
            f2eExpIsDirMat=self.TRI3_F2E_DIR;
            [eMat,f2eMat,f2eIsDirMat] = mapface2edge(vMat,fMat);
            self.assertTrue(isequal(eMat,eExpMat));
            self.assertTrue(isequal(f2eMat,f2eExpMat));
            self.assertTrue(isequal(f2eIsDirMat,f2eExpIsDirMat));
        end
        %
        function testSphereTriExt(self)
            N_REQUESTED_POINTS = 500;
            N_FACES_EXPECTED=1280;
            N_VERTICES_EXPECTED=642;
            [vMat,fMat] = spheretri(N_REQUESTED_POINTS);
            self.assertEqual(size(vMat,1),N_VERTICES_EXPECTED);
            self.assertEqual(size(fMat,1),N_FACES_EXPECTED);
            try
                patch('Vertices',vMat,'Faces',fMat,'FaceColor','g',...
                    'EdgeColor','black');
                close all;
            catch meObj
                close all;
                rethrow(meObj);
            end
        end
        function testShrinkFaceTri(self)
            MAX_DIST=0.5;
            N_DATA_SETS=2;
            for iDataSet=N_DATA_SETS:-1:1
                SInp=self.loadData(['inp',num2str(iDataSet)]);
                [v0Mat,f0Mat]=deal(SInp.v0,SInp.f0);
                [v1Mat,f1Mat]=shrink(v0Mat,f0Mat);
                %% check that no vertices is deleted
                isOldVertsKept=all(ismember(v0Mat,v1Mat,'rows'));
                self.assertTrue(isOldVertsKept);
                %% check that all edges are short enough
                tr=triangulation(f1Mat,v1Mat);
                e1Mat=tr.edges();
                dMat=v1Mat(e1Mat(:,1),:)-v1Mat(e1Mat(:,2),:);
                maxEdgeLength=max(realsqrt(sum(dMat.*dMat,2)));
                self.assertTrue(maxEdgeLength<=MAX_DIST);
                %% regression test
                [SOut.v0,SOut.f0]=deal(v1Mat,f1Mat);
                %
                %self.saveData(['out',num2str(iDataSet)],SOut);
                %
                SEOut=self.loadData(['out',num2str(iDataSet)]);
                self.assertTrue(isequal(SOut,SEOut));
            end
            function [v1Mat,f1Mat]=shrink(v0Mat,f0Mat)
                import gras.geom.tri.*;
                %% shrink faces
                [v1Mat,f1Mat,S1Stat]=self.aux_shrinkfacetri(v0Mat,...
                    f0Mat,MAX_DIST);
                %% Perform additional checks
                [v2Mat,f2Mat,S2Stat]=self.aux_shrinkfacetri(v0Mat,...
                    f0Mat,MAX_DIST,S1Stat.nSteps);
                self.assertTrue(isequal(v1Mat,v2Mat));
                self.assertTrue(isequal(f1Mat,f2Mat));
                self.assertTrue(isequal(S1Stat,S2Stat));
                %
                checkStepWise(v0Mat,f0Mat,0,3);
                checkStepWise(v0Mat,f0Mat,MAX_DIST);
            end
            function checkStepWise(vInpMat,fInpMat,maxTol,varargin)
                MAX_TOL=1e-14;
                [vResMat,fResMat,SStat]=self.aux_shrinkfacetri(vInpMat,...
                    fInpMat,maxTol,varargin{:});
                nSteps=SStat.nSteps;
                if nSteps>1
                    [v1Mat,f1Mat,~]=self.aux_shrinkfacetri(vInpMat,...
                        fInpMat,maxTol,nSteps-1);
                    [v2Mat,f2Mat,~]=self.aux_shrinkfacetri(v1Mat,...
                        f1Mat,maxTol,1);
                    [isPos,reportStr]=istriequal(vResMat,fResMat,...
                        v2Mat,f2Mat,MAX_TOL);
                    self.assertTrue(isPos,reportStr);
                end
            end
        end
        %
    end
end