;;----- This procedure borrowed heavily from CWO's /home/codell/idlprogs/aos/models/OCO/l2/b9_acos/plot_acos_vs_models_b9.pro
;;----- Modifications by T.E.T Autumn, 2020 to adapt and generate analysis/plots for the GOSAT v9 data paper.
;;----- New version <main_satellite_vs_model.pro> developed for the OCO ROSES Tropical IAV work to harmonize GOSAT v9 and OCO-2 v10.
;;-----
;;----- Keep in mind that the matchup of models to OCO uses the ten second averages. So the "sounding densities" will be 
;;----- roughly order 10sec*3Hz*8sound/frame=240 times smaller than native. So instead of approx 100k good quality flagged soundings per day
;;----- there will be more like 400 matchups per day. Translates to 150k matchups per year.


;;----- Set the sensor and version
IF (1) THEN BEGIN
   VerName_1='gosat_v9a'
   VerName_2='gosat_v9c2'
   mode_list=['landH','landM','oceanH']
   ;mode_list=['landH']
   ;mode_list=['oceanH']
ENDIF

Col_1a='dark green'
Col_1b='lime'
Col_2a='orange'
Col_2b='red'
Col_3a='blue'
Col_3b='aqua'
Title='GOSAT - MMM (binned at 30 days)'

;;----- Plot settings.
Margin=[0.1,0.2,0.1,0.2]
NCols=1 & NRows=3
FontSize_1=10
FontSize_2=8
FontStyle=1
FitThick=2
BinSymSize=0.75
Symbol='Circle'
SymIncrement=1


FOR mode=0,2 DO BEGIN

      CASE Mode OF
         0: BEGIN
               PlotPos=1
               CURRENT=0
               PlotLabel='(a)'
               XTitle=''
               YTitle=''
               Color_1=Col_1a
               Color_2=Col_1b
            END
         1: BEGIN
               PlotPos=2
               CURRENT=1
               PlotLabel='(b)'
               XTitle=''
               YTitle='$\Delta$XCO$_2$ [ppm] (GOSAT $-$ MMM)'
               Color_1=Col_2a
               Color_2=Col_2b
            END
         2: BEGIN
               PlotPos=3
               CURRENT=1
               PlotLabel='(c)'
               XTitle=''
               YTitle=''
               Color_1=Col_3a
               Color_2=Col_3b
            END
      ENDCASE

   InFile_1='/home/ttaylor/analysis_utilities/tropical_iav/plots/satellite_vs_models/gosat_v9a/gosat_v9a_vs_4models_bc_ak_binned_Global_'+mode_list[mode]+'_20090420-20200630_Annual.sav'
   InFile_2='/home/ttaylor/analysis_utilities/tropical_iav/plots/satellite_vs_models/gosat_v9c2/gosat_v9c2_vs_4models_bc_ak_binned_Global_'+mode_list[mode]+'_20090420-20200630_Annual.sav'

   RESTORE,FILENAME=InFile_1,VERBOSE=1
   dxco2_bin_v9a=dxco2_bin
   RESTORE,FILENAME=InFile_2,VERBOSE=1
   dxco2_bin_v9c=dxco2_bin

   YRange=[-1.0,1.0]

   ;;----- X axis setup.
   Dummy = LABEL_DATE(DATE_FORMAT='%N/%Z')
   XTitle='Time [MM/YY]'
   stretch_lower= ( MAX(dxco2_bin_v9a.xb) - MIN(dxco2_bin_v9a.xb) ) * 0.05
   stretch_upper= ( MAX(dxco2_bin_v9a.xb) - MIN(dxco2_bin_v9a.xb) ) * 0.25
   xrange_plot=[MIN(dxco2_bin_v9a.xb)-stretch_lower,MAX(dxco2_bin_v9a.xb)+stretch_upper]

   ;;----- Create initial plot
   p1=PLOT(dxco2_bin_v9a.xb, dxco2_bin_v9a.yb, $
           Symbol='o',Sym_Size=1,Sym_Filled=1,Sym_color=Color_1,Sym_Fill_Color=Color_1,Color=Color_1, $
           Title=Title,XTitle='Julian Day',YTitle='$\Delta$XCO$_2$ (ppm)', $
           YRange=YRange, XRange=xrange_plot, XTickFormat='LABEL_DATE', $
           NAME=VerName_1, CURRENT=Current, LAYOUT=[NCols,NRows,PlotPos], MARGIN=Margin $
)

   p2=PLOT(OVERPLOT=p1,dxco2_bin_v9c.xb,dxco2_bin_v9c.yb,Symbol='d',Sym_Size=1,Sym_Filled=1,Sym_color=Color_2,Sym_Fill_Color=Color_2,Color=Color_2,NAME=VerName_2)

   ;;----- .compile /usr/local/harris/idl/lib/graphics/legend.pro
   text = TEXT(0.85,0.80,VerName_1,TARGET=p1,/RELATIVE, FONT_SIZE=FontSize_2,FONT_COLOR=Color_1,FONT_STYLE=FontStyle)
   text = TEXT(0.85,0.65,VerName_2,TARGET=p2, /RELATIVE, FONT_SIZE=FontSize_2,FONT_COLOR=Color_2,FONT_STYLE=FontStyle)

   ;;----- Overplot horizontal zero lines
   p3=PLOT(OVERPLOT=p2,XRange_plot,[0.,0.],Thick=2,Color=Black)

ENDFOR

STOP          
END
