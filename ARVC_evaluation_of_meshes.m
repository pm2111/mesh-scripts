volumes = [];
distances = {};
subjects = {'ARVC_002','ARVC_005','ARVC_009','ARVC_018'};
cd D:\ARVC_MESHES\
for i=1:4
subject = subjects{i};

HQ = loadv( fullfile(subject,'HS'),'HS');
HQ = contoursFrom(  HQ , loadv( fullfile(subject,'HQ'),'HQ') );
HQ = transformFrom( HQ , loadv( fullfile(subject,'HQ'),'HQ') );


iZ = loadv( fullfile(subject,'HM'),'iZ');
HQ = transform( HQ , iZ );

EPIm = transform( loadv( fullfile(subject,'HM') , 'EPIms' ) , iZ );
LVm  = transform( loadv( fullfile(subject,'HM') , 'LVms'  ) , iZ );
RVm  = transform( loadv( fullfile(subject,'HM') , 'RVms'  ) , iZ );


%%

M = EPIm; C = {};
for r = 1:size( HQ ,1 )
  tC = [ HQ{r,2} ; HQ{r,5} ];
  if ~isempty(tC)
    tC = double( resample( polyline( tC ) ,'e',0.5) );
    w = tC(:,3) > 0; tC(w,:) = [];
  end
  C{r,1} = tC( ~any(isnan(tC),2) ,:);
end
C = cell2mat( C );
[~,cp,d] = vtkClosestElement( Mesh(M,0) , C );
volumes(i,1)=meshVolume( MeshFillHoles( LVm ,Inf) );
volumes(i,2)=meshVolume( MeshFillHoles( EPIm ,Inf) );

subplot(1,2,1)
plotMESH( M , 'FaceAlpha',0.4,'ne','gouraud');
hplot3d( C , cp , '-','color',[1,1,1]*.5)
hold on; cline( C , d , 'Marker','o','MarkerSize',4,'MarkerFaceColor','flat','MarkerEdgeColor','none'); hold off;
axis(objbounds(gca,1.1));
colormap jet
caxis( [ 0 , prctile( d ,90 ) ] );

subplot(1,2,2)
hist(d,50); vline(   prctile( d ,90 ));
distances{i,1} = d;
%%
M = LVm; C = {};
for r = 1:size( HQ ,1 )
  tC = [ HQ{r,3} ];
  if ~isempty(tC)
    tC = double( resample( polyline( tC ) ,'e',0.5) );
    w = tC(:,3) > 0; tC(w,:) = [];
  end
  C{r,1} = tC( ~any(isnan(tC),2) ,:);
end
C = cell2mat( C );
[~,cp,d] = vtkClosestElement( Mesh(M,0) , C );
distances{i,2} = d;

subplot(121)
plotMESH( M , 'FaceAlpha',0.4,'ne','gouraud');
hplot3d( C , cp , '-','color',[1,1,1]*.5)
hold on; cline( C , d , 'Marker','o','MarkerSize',4,'MarkerFaceColor','flat','MarkerEdgeColor','none'); hold off;
axis(objbounds(gca,1.1));
colormap jet
caxis( [ 0 , prctile( d ,90 ) ] );

subplot(122)
hist(d,50); vline(prctile( d ,90 ));

%%

M = RVm; C = {};
for r = 1:size( HQ ,1 )
  tC = [ HQ{r,4} ; HQ{r,6} ];
  if ~isempty(tC)
    tC = double( resample( polyline( tC ) ,'e',0.5) );
    w = tC(:,3) > 0; tC(w,:) = [];
  end
  C{r,1} = tC( ~any(isnan(tC),2) ,:);
end
C = cell2mat( C );
[~,cp,d] = vtkClosestElement( Mesh(M,0) , C );
distances{i,3} = d;

subplot(121)
figure()
plotMESH( M , 'FaceAlpha',0.4,'ne','gouraud');
hplot3d( C , cp , '-','color',[1,1,1]*.5)
hold on; cline( C , d , 'Marker','o','MarkerSize',4,'MarkerFaceColor','flat','MarkerEdgeColor','none'); hold off;
axis(objbounds(gca,1.1));
colormap jet
caxis( [ 0 , prctile( d ,90 ) ] );
colorbar
axis off
title('Point of contour to nearest surface point distance [mm]','Fontsize',35)
set(findobj(gca,'Type','text'),'FontSize',35)
set(gca,'fontSize',35)


subplot(122)
hist(d,50); vline(prctile( d ,90 ));

%%

mean(d),prctile(d,[50,75,90,100])
volumes(i,3)=meshVolume( MeshFillHoles( RVm ,Inf) );
volumes(i,4)=sum( meshQuality( RVm ,'area' ) ) * mean( d ); %error in volume RV

AS = [];
for r = 1:size(HQ,1)
  if ~strcmp( HQ{r,1}.INFO.PlaneName , 'SAx' ), continue; end
  C = [ HQ{r,4} ; HQ{r,6} ]; if isempty( C ), continue; end
  AS(end+1,1) = area( polygon(C(:,1:2)) );
  AS( end ,2) = median( C(:,3) );
end

sum(AS(:,1))*mean( diff( AS(:,2) ) )
end