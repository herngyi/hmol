module Main where

import Geometry
import Chemistry
import MoleculeBuilder
import Rendering
import Molecules

import Algorithms
import Proteins
import NucleicAcids
import NanoPutians

import Graphics.UI.Gtk hiding (Point, Fill, get, Color)
import qualified Graphics.UI.Gtk as Gtk
import Graphics.UI.Gtk.Gdk.EventM
import Graphics.UI.Gtk.Glade
import Graphics.UI.Gtk.OpenGL
import Graphics.Rendering.OpenGL

import Data.List
import Data.Array
import Data.IORef
import Control.Monad
import Control.Monad.Trans
import Control.Exception
import System.Directory

tr vm = moleculeGetRends vm (nullToggleOptions // [(showConvexHulls, True)]) Allen

data Project = Project String (Maybe FilePath) Molecule RendMolecule ViewMode ToggleOptions ElectronegativityType Double

nullProject :: String -> Project
nullProject str = Project str Nothing nullMolecule nullRendMol Default nullToggleOptions Allen 1

main :: IO ()
main = do
  unsafeInitGUIForThreadedRTS
  initGL
  glconfig <- glConfigNew [GLModeRGBA, GLModeDepth, GLModeDouble]

  Just xml <- xmlNew "Interface.glade"

  window <- xmlGetWidget xml castToWindow "Main"
  new    <- xmlGetWidget xml castToImageMenuItem "imagemenuitem3"
  open   <- xmlGetWidget xml castToImageMenuItem "imagemenuitem1"
  save   <- xmlGetWidget xml castToImageMenuItem "imagemenuitem2"
  saveAs <- xmlGetWidget xml castToImageMenuItem "imagemenuitem4"
  quit   <- xmlGetWidget xml castToImageMenuItem "imagemenuitem5"
  onDestroy window mainQuit
  onActivateLeaf quit (widgetDestroy window)

  openButton  <- xmlGetWidget xml castToButton "button9"
  saveButton  <- xmlGetWidget xml castToButton "button10"
  closeButton <- xmlGetWidget xml castToButton "button8"
  widgetSetSensitivity saveButton  False
  widgetSetSensitivity closeButton False

  filePath <- newIORef Nothing

  erWindow <- xmlGetWidget xml castToWindow "Error"
  erLabel  <- xmlGetWidget xml castToLabel  "label3"
  erOk     <- xmlGetWidget xml castToButton "button5"
  onClicked erOk (widgetHide erWindow)
  let errorMessage title message = do
        (win:_) <- filterM windowHasToplevelFocus =<< windowListToplevels
        windowSetTransientFor erWindow win
        windowSetTitle erWindow title
        labelSetText erLabel message
        widgetShow erWindow

  cfWindow <- xmlGetWidget xml castToWindow "Confirmation"
  cfLabel  <- xmlGetWidget xml castToLabel  "label15"
  cfYes    <- xmlGetWidget xml castToButton "button6"
  cfNo     <- xmlGetWidget xml castToButton "button7"
  cfAction <- newIORef (return ())
  onClicked cfYes $ do
    action <- readIORef cfAction
    widgetHide cfWindow
    action
  onClicked cfNo (widgetHide cfWindow)
  let confirm title message action = do
        writeIORef cfAction action
        win:_ <- filterM windowHasToplevelFocus =<< windowListToplevels
        windowSetTransientFor cfWindow win
        windowSetTitle cfWindow title
        labelSetText cfLabel message
        widgetShow cfWindow

  currMol  <- newIORef nullMolecule
  rendMol  <- newIORef nullRendMol
  viewMode <- newIORef Default
  toggOpts <- newIORef nullToggleOptions
  enType   <- newIORef Allen

  translation <- newIORef origin
  rotation    <- newIORef identity
  scaleFactor <- newIORef 1

  canvas <- glDrawingAreaNew glconfig
  widgetAddEvents canvas [PointerMotionMask, ScrollMask, ButtonPressMask]

  canvasPort <- xmlGetWidget xml castToViewport "viewport1"
  containerAdd canvasPort canvas

  onRealize canvas $ withGLDrawingArea canvas $ \_ -> do
    clearColor $= Color4 1 1 1 1
    perspective 0.785 1 (-1) 1
    lookAt (Vertex3 0 0 (-1)) (Vertex3 0 0 0) (Vector3 0 1 0)
    depthFunc  $= Just Less
    drawBuffer $= BackBuffers

    attenuation (Light 0) $= (1,0,0)

    position (Light 0) $= Vertex4 1 1 1 0
    lightModelAmbient $= Color4 0.5 0.5 0.5 1
    shadeModel $= Smooth

    normalize $= Enabled

    ambient (Light 0) $= Color4 0.2 0.2 0.2 1
    diffuse (Light 0) $= Color4 0.5 0.5 0.5 1
    specular (Light 0) $= Color4 0.5 0.5 0.5 1

    colorMaterial $= Just (Front, AmbientAndDiffuse)

--    materialSpecular Front $= Color4 1 1 1 1
--    materialShininess Front $= 128

    light (Light 0) $= Enabled
    lighting        $= Enabled

    blendFunc $= (SrcAlpha, OneMinusSrcAlpha)
    blend $= Enabled

  orientation  <- newIORef 0
  orientScroll <- xmlGetWidget xml castToVScrollbar "vscrollbar1"

  onRangeValueChanged orientScroll $ do
    t <- rangeGetValue orientScroll
    s <- readIORef orientation
    writeIORef orientation t
    modifyIORef rendMol (rotate3D (Vector 0 0 1) (s - t))
    widgetQueueDraw canvas

  let recalculate = do
        mol <- readIORef currMol
        vm  <- readIORef viewMode
        tgs <- readIORef toggOpts
        ent <- readIORef enType
        let rends@(RendMolecule size _ _ _) = moleculeGetRends vm tgs ent mol
        writeIORef rendMol rends
        (w,h) <- widgetGetSize canvas
        let k = if w < h
                 then fromIntegral w / fromIntegral h
                 else fromIntegral h / fromIntegral w
        writeIORef scaleFactor (k / size)
        widgetQueueDraw canvas
        writeIORef orientation 0
        rangeSetValue orientScroll 0

  mass     <- xmlGetWidget xml castToLabel "label6"
  charge   <- xmlGetWidget xml castToLabel "label10"
  diameter <- xmlGetWidget xml castToLabel "label12"

  crystalData    <- xmlGetWidget xml castToVBox  "vbox13"
  crystalDensity <- xmlGetWidget xml castToLabel "label18"
  crystalA       <- xmlGetWidget xml castToLabel "label20"
  crystalB       <- xmlGetWidget xml castToLabel "label25"
  crystalC       <- xmlGetWidget xml castToLabel "label26"
  crystalAlpha   <- xmlGetWidget xml castToLabel "label30"
  crystalBeta    <- xmlGetWidget xml castToLabel "label31"
  crystalGamma   <- xmlGetWidget xml castToLabel "label32"

  let renewSpecialData (CrystalData n dsy [ja,jb,jc,jalpha,jbeta,jgamma]) = do
        labelSetText crystalDensity (show (roundMass dsy) ++ suffix n)
        labelSetText crystalA     (maybe "= -" showLength ja    )
        labelSetText crystalB     (maybe "= -" showLength jb    )
        labelSetText crystalC     (maybe "= -" showLength jc    )
        labelSetText crystalAlpha (maybe "= -" showAngle  jalpha)
        labelSetText crystalBeta  (maybe "= -" showAngle  jbeta )
        labelSetText crystalGamma (maybe "= -" showAngle  jgamma)
        widgetShow crystalData
        where suffix 1 = " g/cm"
              suffix n = " g/cm^" ++ show n
              showLength p = "= " ++ show (roundDiameter p) ++ " pm"
              showAngle  p = "= " ++ show (roundAngle    p) ++ " rad"
              roundAngle x = fromIntegral (round (x * 10^3)) / 10^3

      renewData mol = do
        let prs = map (\(Atom p e _ _ _) -> (p, atomicRadii ! e)) (vertices mol)
            MolData m q sds = molData mol
            d = maximum (0 : distances prs)
              where distances (pr : prs) = map (getDists pr) prs ++ distances prs
                    distances _ = []
                    getDists (p1,r1) (p2,r2) = distance p1 p2 + r1 + r2
        labelSetText mass     (show (roundMass m) ++ " g/mol")
        labelSetText charge   (showCharge q)
        labelSetText diameter (show (roundDiameter d) ++ " pm")
        widgetHide crystalData
        mapM_ renewSpecialData sds

  vmDefault <- xmlGetWidget xml castToRadioButton "radiobutton1"
  onToggled vmDefault $ do
    active <- toggleButtonGetActive vmDefault
    when active (writeIORef viewMode  Default >> recalculate)

  vmSpaceFilling <- xmlGetWidget xml castToRadioButton "radiobutton2"
  onToggled vmSpaceFilling $ do
    active <- toggleButtonGetActive vmSpaceFilling
    when active (writeIORef viewMode  SpaceFilling >> recalculate)

  vmSkeletal <- xmlGetWidget xml castToRadioButton "radiobutton3"
  onToggled vmSkeletal $ do
    active <- toggleButtonGetActive vmSkeletal
    when active (writeIORef viewMode  Skeletal >> recalculate)

  opContinuousBondColoring <- xmlGetWidget xml castToCheckButton "checkbutton4"
  onToggled opContinuousBondColoring $ do
    active <- toggleButtonGetActive opContinuousBondColoring
    modifyIORef toggOpts (// [(continuousBondColoring, active)])
    recalculate

  opShowLonePairs <- xmlGetWidget xml castToCheckButton "checkbutton3"
  onToggled opShowLonePairs $ do
    active <- toggleButtonGetActive opShowLonePairs
    modifyIORef toggOpts (// [(showLonePairs, active)])
    recalculate

  opShowEmptyOrbitals <- xmlGetWidget xml castToCheckButton "checkbutton6"
  onToggled opShowEmptyOrbitals $ do
    active <- toggleButtonGetActive opShowEmptyOrbitals
    modifyIORef toggOpts (// [(showEmptyOrbitals, active)])
    recalculate

  opShowConvexHulls <- xmlGetWidget xml castToCheckButton "checkbutton8"
  onToggled opShowConvexHulls $ do
    active <- toggleButtonGetActive opShowConvexHulls
    modifyIORef toggOpts (// [(showConvexHulls, active)])
    recalculate

  flMassDensity <- xmlGetWidget xml castToCheckButton "checkbutton1"
  onToggled flMassDensity $ do
    active <- toggleButtonGetActive flMassDensity
    modifyIORef toggOpts (// [(filterMassDensity, active)])
    recalculate

  flChargeDensity <- xmlGetWidget xml castToCheckButton "checkbutton7"
  onToggled flChargeDensity $ do
    active <- toggleButtonGetActive flChargeDensity
    modifyIORef toggOpts (// [(filterChargeDensity, active)])
    recalculate

  flElectronegativity <- xmlGetWidget xml castToCheckButton "checkbutton5"
  onToggled flElectronegativity $ do
    active <- toggleButtonGetActive flElectronegativity
    modifyIORef toggOpts (// [(filterElectronegativity, active)])
    recalculate

  enCombo <- xmlGetWidget xml castToComboBox "combobox1"
  comboBoxSetActive enCombo 0
  on enCombo changed $ do
    n <- comboBoxGetActive enCombo
    let ent = case n of
                0 -> Allen
                1 -> Pauling
                2 -> AllredRochow
                3 -> Mulliken
                4 -> Jaffe
                5 -> Sanderson
    writeIORef enType ent
    recalculate

  dragButton    <- newIORef Nothing
  dragTranslate <- newIORef Nothing
  dragRotate    <- newIORef Nothing

  onExpose canvas $ \_ -> do
    withGLDrawingArea canvas $ \glwindow -> do
      (w,h) <- widgetGetSize canvas
      let s = fromIntegral (max w h)
          (x,y) = if w < h
                   then ((fromIntegral w - fromIntegral s) / 2, 0)
                   else (0, (fromIntegral h - fromIntegral s) / 2)
      viewport $= (Position (round x) (round y), Size s s)
      clear [DepthBuffer, ColorBuffer]
      rm <- readIORef rendMol
      trans <- readIORef translation
      (axis, ang) <- readIORef rotation
      scl <- readIORef scaleFactor
      let io = do
            scaleGL scl
            rotateRadGL ang axis
            translateGL trans
          a  = 1 / scl
          v1 = translate3D (-trans) (rotate3D axis (-ang) (Vector a a a))
          b  = -a
          v2 = translate3D (-trans) (rotate3D axis (-ang) (Vector b b b))
          ns = map (rotate3D axis (-ang)) [Vector 1 0 0, Vector 0 1 0, Vector 0 0 1]
          xyz = map (\n -> sort [n `dot` v1, n `dot` v2]) ns
      renderMolecule io (fromIntegral s * scl) ns xyz rm
      glDrawableSwapBuffers glwindow
    return True

  on canvas buttonPressEvent $ do
    (x,y) <- eventCoordinates
    bt    <- eventButton
    liftIO $ case bt of
               RightButton -> do
                 writeIORef dragButton (Just RightButton)
                 writeIORef dragTranslate $ Just (x,y)
               LeftButton -> do
                 writeIORef dragButton (Just LeftButton)
                 (w,h) <- widgetGetSize canvas
                 let k = pi / fromIntegral (max w h)
                 writeIORef dragRotate $ Just ((x,y), k)
    return True

  on canvas buttonReleaseEvent $ do
    bt <- eventButton
    liftIO $ do
      case bt of
        RightButton -> do
          trans <- readIORef translation
          modifyIORef rendMol (translate3D trans)
          writeIORef dragTranslate Nothing
          writeIORef translation origin
        LeftButton -> do
          (axis, ang) <- readIORef rotation
          modifyIORef rendMol (rotate3D axis ang)
          writeIORef dragRotate Nothing
          writeIORef rotation identity
      writeIORef dragButton Nothing
      widgetQueueDraw canvas
    return True

  on canvas motionNotifyEvent $ do
    (px,py) <- eventCoordinates
    liftIO $ do
      bt <- readIORef dragButton
      case bt of
        Just RightButton -> liftIO $ do
          (w,h) <- widgetGetSize canvas
          scl <- readIORef scaleFactor
          Just (ox,oy) <- readIORef dragTranslate
          let k = scl * fromIntegral (max w h)
              v = (2 / k) |*| Vector (px - ox) (oy - py) 0
          writeIORef translation v
          widgetQueueDraw canvas
          return True
        Just LeftButton -> liftIO $ do
          Just ((ox,oy), k) <- readIORef dragRotate
          let ndx  = ox - px
              ndy  = oy - py
              l    = sqrt (ndx^2 + ndy^2)
              axis = Vector (ndy / l) (ndx / l) 0
              ang  = k * l
          writeIORef rotation (axis, ang)
          widgetQueueDraw canvas
          return True
        _ -> return False

  on canvas scrollEvent $ do
    dir <- eventScrollDirection
    liftIO $ do
      scl <- readIORef scaleFactor
      let scl' = case dir of
                   ScrollUp   -> scl * 1.1
                   ScrollDown -> scl / 1.1
      writeIORef scaleFactor scl'
      widgetQueueDraw canvas
    return True

  projects <- newIORef []
  untitled <- newIORef 1
  currProj <- newIORef (-1)

  projToolbar <- xmlGetWidget xml castToToolbar "toolbar1"
  
  invisible   <- radioToolButtonNew
  projToggles <- newIORef []

  let switchToProject n = do
        m <- readIORef currProj
        unless (m == n) $ do
          pjs <- readIORef projects
          let Project str _ _ _ _ _ _ _ = pjs !! m
          jfp <- readIORef filePath
          mol <- readIORef currMol
          rs  <- readIORef rendMol
          vm  <- readIORef viewMode
          tgs <- readIORef toggOpts
          ent <- readIORef enType
          scl <- readIORef scaleFactor
          let (as,_:bs) = splitAt m pjs
          writeIORef projects (as ++ (Project str jfp mol rs vm tgs ent scl : bs))        

          let Project str' jfp' mol' rs' vm' tgs' ent' scl' = pjs !! n
          writeIORef currMol nullMolecule
          setViewOptions vm' tgs' ent'
          renewData  mol'
          writeIORef filePath jfp'
          writeIORef currMol mol'
          writeIORef rendMol rs'
          writeIORef viewMode vm'
          writeIORef toggOpts tgs'
          writeIORef enType ent'
          writeIORef scaleFactor scl'
          writeIORef currProj n
          widgetQueueDraw canvas
          widgetSetSensitivity saveButton  True
          widgetSetSensitivity closeButton True

      setViewOptions vm tgs ent = do
        flip toggleButtonSetActive True $
          case vm of
            Default      -> vmDefault
            SpaceFilling -> vmSpaceFilling
            Skeletal     -> vmSkeletal
        toggleButtonSetActive flMassDensity            (tgs ! filterMassDensity)
        toggleButtonSetActive flChargeDensity          (tgs ! filterChargeDensity)
        toggleButtonSetActive flElectronegativity      (tgs ! filterElectronegativity)
        toggleButtonSetActive opContinuousBondColoring (tgs ! continuousBondColoring)
        toggleButtonSetActive opShowLonePairs          (tgs ! showLonePairs)
        toggleButtonSetActive opShowEmptyOrbitals      (tgs ! showEmptyOrbitals)
        toggleButtonSetActive opShowConvexHulls        (tgs ! showConvexHulls)

      newProject str = do
        tgs <- readIORef projToggles
        let n = length tgs
        newRadio <- radioToolButtonNewFromWidget invisible
        let newToggle = toToggleToolButton newRadio
        toolButtonSetLabel newToggle (Just str)
        toolbarInsert projToolbar newToggle (-1)
        set projToolbar [toolbarChildHomogeneous newToggle := False]
        widgetShow newToggle
        onToolButtonToggled newToggle (switchToProject n)
        tgs <- readIORef projToggles
        writeIORef projToggles (tgs ++ [Just newToggle])
        modifyIORef projects (++ [nullProject str])
        toggleToolButtonSetActive newToggle True
        widgetSetSensitivity saveButton  True
        widgetSetSensitivity closeButton True

  onActivateLeaf new $ do
    n <- readIORef untitled
    newProject ("Molecule " ++ show n)
    writeIORef untitled (n + 1)

  openChooser <- xmlGetWidget xml castToWindow "window1"
  opOk     <- xmlGetWidget xml castToButton "button3"
  opCancel <- xmlGetWidget xml castToButton "button4"
  onClicked opCancel (widgetHide openChooser)

  openVBox <- xmlGetWidget xml castToVBox "vbox8"
  openMol  <- fileChooserWidgetNew FileChooserActionOpen
  boxPackStart openVBox openMol PackGrow 0
  opFilter <- fileFilterNew
  fileFilterSetName opFilter "Molecule Data File (*.mol)"
  fileFilterAddPattern opFilter "*.mol"
  fileChooserAddFilter openMol opFilter

  onActivateLeaf open  (widgetShowAll openChooser)
  onClicked openButton (widgetShowAll openChooser)

  afterFileActivated openMol $ do
    jfn@(Just filename) <- fileChooserGetFilename openMol
    str  <- readFile filename
    emol <- try (return $! read str) :: IO (Either ErrorCall Molecule)
    case emol of
      Right mol -> do
        widgetHide openChooser
        let short = takeWhile (/= '.') (reverse (takeWhile (/= '\\') (reverse filename)))
        newProject short
        writeIORef currMol mol
        recalculate
        renewData mol
        writeIORef filePath jfn
      Left e -> do
        let fn = reverse (takeWhile (/= '\\') (reverse filename))
        errorMessage "Molecule Parsing Error" ('\"' : fn ++ "\" could not be parsed into a molecule.")

  onClicked opOk $ do
    jfn <- fileChooserGetFilename openMol
    case jfn of
      Just filename -> do
        str <- readFile filename
        emol <- try (return $! read str) :: IO (Either ErrorCall Molecule)
        case emol of
          Right mol -> do
            widgetHide openChooser
            let short = takeWhile (/= '.') (reverse (takeWhile (/= '\\') (reverse filename)))
            newProject short
            writeIORef currMol mol
            recalculate
            renewData mol
            writeIORef filePath jfn
          Left e -> do
            let fn = reverse (takeWhile (/= '\\') (reverse filename))
            errorMessage "Molecule Parsing Error" ('\"' : fn ++ "\" could not be parsed into a molecule.")
      _ -> return ()

  saveChooser <- xmlGetWidget xml castToWindow "window2"
  svOk     <- xmlGetWidget xml castToButton "button1"
  svCancel <- xmlGetWidget xml castToButton "button2"
  onClicked svCancel (widgetHide saveChooser)

  saveVBox <- xmlGetWidget xml castToVBox "vbox3"
  saveMol  <- fileChooserWidgetNew FileChooserActionSave
  boxPackStart saveVBox saveMol PackGrow 0
  fileChooserSetDoOverwriteConfirmation saveMol True

  onActivateLeaf saveAs (widgetShowAll saveChooser)
  onActivateLeaf save $ do
    jfp <- readIORef filePath
    case jfp of
      Just fp -> do
        mol <- readIORef currMol
        writeFile fp (show mol)
      _ -> widgetShowAll saveChooser
  onClicked saveButton $ do
    jfp <- readIORef filePath
    case jfp of
      Just fp -> do
        mol <- readIORef currMol
        writeFile fp (show mol)
      _ -> widgetShowAll saveChooser

  onClicked svOk $ do
    jfns <- fileChooserGetFilename saveMol
    case jfns of
      Just fn -> do
        mol <- readIORef currMol
        let filename = takeWhile (/= '.') fn ++ ".mol"
            short = reverse (takeWhile (/= '\\') (reverse filename))
        overwrite <- doesFileExist filename
        if overwrite
         then confirm "Overwrite Existing File" ('\"' : short ++ "\" already exists. Are you sure you want to overwrite it?") $ do
                n <- readIORef currProj
                tgs <- readIORef projToggles
                let Just tg = tgs !! n
                toolButtonSetLabel tg (Just (takeWhile (/= '.') short))
                writeFile  filename (show mol)
                writeIORef filePath (Just filename)
                widgetHide saveChooser
         else do
           n <- readIORef currProj
           tgs <- readIORef projToggles
           let Just tg = tgs !! n
           toolButtonSetLabel tg (Just (takeWhile (/= '.') short))
           writeFile  filename (show mol)
           writeIORef filePath (Just filename)
           widgetHide saveChooser
      _ -> return ()

  let closeProject = do
        n   <- readIORef currProj
        tgs <- readIORef projToggles
        let (as,(Just tg):bs) = splitAt n tgs
        widgetHide tg
        writeIORef projToggles (as ++ (Nothing : bs))
        case find (/= Nothing) (as ++ bs) of
          Just (Just tg') -> toggleToolButtonSetActive tg' True
          _ -> do
            writeIORef filePath Nothing
            writeIORef currMol nullMolecule
            writeIORef rendMol nullRendMol
            writeIORef viewMode Default
            writeIORef toggOpts nullToggleOptions
            writeIORef enType Allen
            writeIORef scaleFactor 1
            writeIORef currProj  (-1)
            widgetQueueDraw canvas
            widgetSetSensitivity saveButton  False
            widgetSetSensitivity closeButton False

  onClicked closeButton closeProject

  on window keyPressEvent $ do
    key <- eventKeyName
    liftIO $ do
      case key of
        "Shift_L" -> return True
        _         -> return False

  on window keyReleaseEvent $ do
    key <- eventKeyName
    liftIO $ do
      case key of
        "Shift_L" -> return True
        _         -> return False

  writeIORef currMol nanoBalletDancer
  recalculate
  windowMaximize window
  widgetShowAll window
  mainGUI

t mol folder name = writeFile ("Molecules/" ++ folder ++ "/" ++ name ++ ".mol") (show mol)

full = do
  regenA
  regenB
  regenProt
  regenC
  regenD
  regenE

regenA = do
  t alpha_glucose_ring "Sugars" "alpha-Glucose (Ring)"
  t beta_glucose_ring  "Sugars" "beta-Glucose (Ring)"

  t ethene            "Organic Chemistry" "Ethene"
  t (alkane 10)       "Organic Chemistry" "Decane"
  t ethylenediamine   "Organic Chemistry" "Ethylenediamine"
  t benzene           "Organic Chemistry" "Benzene"
  t benzoic_acid      "Organic Chemistry" "Benzoic Acid"
  t cyclohexane_chair "Organic Chemistry" "Cyclohexane (Chair)"
  t cyclohexane_boat  "Organic Chemistry" "Cyclohexane (Boat)"
  t edta              "Organic Chemistry" "EDTA"

  t hydrogen          "Common Substances" "Hydrogen"
  t water             "Common Substances" "Water"
  t carbon_dioxide    "Common Substances" "Carbon Dioxide"
  t ammonia           "Common Substances" "Ammonia"
  t hydrogen_peroxide "Common Substances" "Hydrogen Peroxide"
  t hydrogen_cyanide  "Common Substances" "Hydrogen Cyanide"

  t nitric_acid     "Common Substances" "Nitric Acid"
  t sulfuric_acid   "Common Substances" "Sulfuric Acid"
  t phosphoric_acid "Common Substances" "Phosphoric Acid"

  t _P4    "Natural Occurences of Elements" "P4"
  t _P4O4  "Natural Occurences of Elements" "P4O4"
  t _P4O10 "Natural Occurences of Elements" "P4O10"

regenB = do
  t (basalAminoAcid L)        "Proteins" "Basal Amino Acid (L)"
  t (basalAminoAcid D)        "Proteins" "Basal Amino Acid (D)"
  t (bondingBasalAminoAcid L) "Proteins" "Bonding Basal Amino Acid (L)"
  t (bondingBasalAminoAcid D) "Proteins" "Bonding Basal Amino Acid (D)"
  t alanine       "Proteins" "L-Alanine"
  t arginine      "Proteins" "L-Arginine"
  t asparagine    "Proteins" "L-Asparagine"
  t aspartic_acid "Proteins" "L-Aspartic Acid"
  t cysteine      "Proteins" "L-Cysteine"
  t glutamic_acid "Proteins" "L-Glutamic Acid"
  t glutamine     "Proteins" "L-Glutamine"
  t glycine       "Proteins" "L-Glycine"
  t histidine     "Proteins" "L-Histidine"
  t isoleucine    "Proteins" "L-Isoleucine"
  t leucine       "Proteins" "L-Leucine"
  t lysine        "Proteins" "L-Lysine"
  t methionine    "Proteins" "L-Methionine"
  t phenylalanine "Proteins" "L-Phenylalanine"
  t proline       "Proteins" "L-Proline"
  t serine        "Proteins" "L-Serine"
  t threonine     "Proteins" "L-Threonine"
  t tryptophan    "Proteins" "L-Tryptophan"
  t tyrosine      "Proteins" "L-Tyrosine"
  t valine        "Proteins" "L-Valine"

  t histidine_zwitterion "Proteins" "L-Histidine Zwitterion"

regenProt = t alphabet "Proteins" "Alphabetic Sequence"

regenC = do
  t adenine  "Nucleic Acids" "Adenine"
  t thymine  "Nucleic Acids" "Thymine"
  t uracil   "Nucleic Acids" "Uracil"
  t guanine  "Nucleic Acids" "Guanine"
  t cytosine "Nucleic Acids" "Cytosine"
  t tata_box "Nucleic Acids" "TATA Box"
  t (telomere 6) "Nucleic Acids" "Human Telomeric Repeat"

regenD = do
  t diamond    "Crystals" "Diamond"
  t cubic_NaCl "Crystals" "Salt (NaCl)"
  t ice_Ih     "Crystals" "Ice (Ih)"
  t ice_Ic     "Crystals" "Ice (Ic)"
  t graphene   "Crystals" "Graphene"
  t graphite   "Crystals" "Graphite"

regenE = do
  t amylose   "Polymers" "Starch (Amylose)"
  t cellulose "Polymers" "Cellulose"