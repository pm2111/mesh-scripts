OUTPUT_FILES = {'mpp/Hmapping.mat'};


%%

H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART.vtk' ) );
[EPI0,LV0,RV0,~,MAT0] = HEARTparts( H0 );

H1 = read_VTK( Fullfile( 'HEART.vtk' ) );
[EPI1,LV1,RV1,~,MAT1] = HEARTparts( H1 );

MAT = MAT1 \ MAT0;

%%
EPI = transform( EPI0 , MAT );
LV  = transform( LV0  , MAT );
RV  = transform( RV0  , MAT );

%% ICP
FPS    = @( M , d ) meshFarthestPointSampling( M , 1 , d , Inf , 1 );
Deform = @( M , uW ) Mesh( uW( M.xyz ) , M );
PRCT   = 0.1;

LS = geospace( 1e-6 , 1e10 , 15 );
DS = linspace( 8 , 5 , numel(LS) ); last_d = NaN;
for it = 1:numel(LS)
  if ~isequal( last_d , DS(it) )
    last_d = DS(it);
    epi1 = FPS( EPI1 ,DS(it));
    lv1  = FPS( LV1  ,DS(it));
    rv1  = FPS( RV1  ,DS(it));
  end
  for itt = 1:50
    fprintf('iteration: %2d.%02d , LAMBDA: %g   , d: %g\n', it , itt , LS(it) , DS(it) );
    if rem( itt , 2 )
      fr_e = FPS( EPI ,DS(it));  [~,to_e] = vtkClosestElement( EPI1 , fr_e );
      fr_l = FPS( LV  ,DS(it));  [~,to_l] = vtkClosestElement( LV1  , fr_l );
      fr_r = FPS( RV  ,DS(it));  [~,to_r] = vtkClosestElement( RV1  , fr_r );
    else
      to_e = epi1;              [~,fr_e] = vtkClosestElement( EPI  , to_e );
      to_l = lv1;               [~,fr_l] = vtkClosestElement( LV   , to_l );
      to_r = rv1;               [~,fr_r] = vtkClosestElement( RV   , to_r );
    end
    fr = [ fr_e ; fr_l ; fr_r ];
    to = [ to_e ; to_l ; to_r ];
    
    to = fr + ( to - fr )*PRCT;
    
    uW = InterpolatingSplines( fr , to , [] , 'r' , 'LAMBDA' , LS(it) );
    EPI = Deform( EPI , uW );
    LV  = Deform( LV  , uW );
    RV  = Deform( RV  , uW );

    if 0 && ~rem( itt , 10 )
      clf;
       plotMESH( H1 ,'ne','FaceAlpha',0.4);
      hplotMESH( EPI ,'b','ne','FaceAlpha',0.5);
      hplotMESH( LV  ,'r','ne','FaceAlpha',0.5);
      hplotMESH( RV  ,'g','ne','FaceAlpha',0.5);
      axis(gca,objbounds(gca,1.1)); headlight; material dull; drawnow
    end
  end
end

%% build the mapping

fr = []; to = [];
[~,id] = FarthestPointSampling( EPI0.xyz  , 1 , 4 ); fr = [ fr ; EPI0.xyz(id,:) ]; [~,x] = vtkClosestElement( EPI1 , EPI.xyz(id,:) ); to = [ to ; x ];
[~,id] = FarthestPointSampling(  LV0.xyz  , 1 , 4 ); fr = [ fr ;  LV0.xyz(id,:) ]; [~,x] = vtkClosestElement(  LV1 ,  LV.xyz(id,:) ); to = [ to ; x ];
[~,id] = FarthestPointSampling(  RV0.xyz  , 1 , 4 ); fr = [ fr ;  RV0.xyz(id,:) ]; [~,x] = vtkClosestElement(  RV1 ,  RV.xyz(id,:) ); to = [ to ; x ];

[u,IS] = InterpolatingSplines( fr , to , [] , 'r' );

%% Saving mapping
Save( 'Hmapping.mat' , 'u' ,'IS','EPI','LV','RV');


%%

hFig = Figure();
 plotMESH( H1  ,'ne' );
hplotMESH( EPI ,'FaceColor','b','ne' );
hplotMESH( LV  ,'FaceColor','r','ne' );
hplotMESH( RV  ,'FaceColor','g','ne' );
axis( gca , objbounds(gca,1.1) );
headlight; material dull; drawnow
Savefig


%% mapping LMK
% 
% LMKS = loadv( fullfile( TORSO_MODEL_DIR , 'LMKS' ) , 'LMKS' );
% 
% LVstimulation = u( LMKS.LVstimulation );
% for it = 1:10, [~,LVstimulation] = vtkClosestElement( LV1 , LVstimulation ); end
% 
% RVstimulation = u( LMKS.RVstimulation );
% for it = 1:10, [~,RVstimulation] = vtkClosestElement( RV1 , RVstimulation ); end
% 
% 
% figure('Name',WHERE_AM_I);
% plotMESH( H1 ,'ne'); camlight; material dull
% hplot3d( LVstimulation ,'ko','markerfacecolor','r'); 
% hplot3d( RVstimulation ,'ko','markerfacecolor','g'); 
% drawnow
% 
% Save( 'LVstimulation.mat'  , 'LVstimulation'  );
% Save( 'RVstimulation.mat'  , 'RVstimulation'  );
% 
% 
% % LV_indexes = [ 254188  242868  264185 275075 ]
% % RV_indexes = [ 60775   35855   61622 ]

%% START mpp_epilogue
%% END mpp_epilogue