module Rendering where

import Geometry
import Chemistry
import MoleculeBuilder
import Graphics.Rendering.OpenGL hiding (color, Plane)
import qualified Graphics.Rendering.OpenGL as GL

import Data.List
import Data.Array
import Data.Maybe
import Control.Monad

-- Interfacing between the Geometric structures

orthoGL :: Double -> Double -> Double -> Double -> Double -> Double -> IO ()
orthoGL l r d t f b = ortho (realToFrac l) (realToFrac r) (realToFrac d) (realToFrac t) (realToFrac f) (realToFrac b)

lookAtGL :: Point -> Point -> Vector -> IO ()
lookAtGL (Vector ex ey ez) (Vector cx cy cz) (Vector ux uy uz) =
  lookAt (Vertex3 (realToFrac ex :: GLdouble)
                  (realToFrac ey :: GLdouble)
                  (realToFrac ez :: GLdouble))
         (Vertex3 (realToFrac cx :: GLdouble)
                  (realToFrac cy :: GLdouble)
                  (realToFrac cz :: GLdouble))
         (Vector3 (realToFrac ux :: GLdouble)
                  (realToFrac uy :: GLdouble)
                  (realToFrac uz :: GLdouble))

vertexGL :: Vector -> IO ()
vertexGL (Vector x y z) = vertex $ Vertex3 (realToFrac x :: GLdouble)
                                           (realToFrac y :: GLdouble)
                                           (realToFrac z :: GLdouble)

normalGL :: Vector -> IO ()
normalGL (Vector x y z) = normal $ Normal3 (realToFrac x :: GLdouble)
                                           (realToFrac y :: GLdouble)
                                           (realToFrac z :: GLdouble)

colorGL :: ColorRGBA -> IO ()
colorGL (ColorRGBA r g b a) = GL.color $ Color4 (realToFrac r :: GLfloat)
                                                (realToFrac g :: GLfloat)
                                                (realToFrac b :: GLfloat)
                                                (realToFrac a :: GLfloat)

translateGL :: Vector -> IO ()
translateGL (Vector x y z) = translate $ Vector3 (realToFrac x :: GLdouble)
                                                 (realToFrac y :: GLdouble)
                                                 (realToFrac z :: GLdouble)
scaleGL :: Double -> IO ()
scaleGL x = scale k k k
  where k = realToFrac x :: GLdouble

rotateRadGL :: Double -> Vector -> IO ()
rotateRadGL rad (Vector x y z) = rotate (realToFrac $ rad * 180 / pi) $
                                          Vector3 (realToFrac x :: GLdouble)
                                                  (realToFrac y :: GLdouble)
                                                  (realToFrac z :: GLdouble)

-- Rendering Data Structures

data ViewMode = Default | SpaceFilling | Skeletal | QuantumOrbitals
  deriving (Show, Eq, Ord, Enum)

defaultScaleRadius :: Double
defaultScaleRadius = 1/6

data Filter = MassDensityFilter | ElectronegativityFilter ElectronegativityType | ElectronDensityFilter
  deriving (Show, Eq)

type ToggleOptions = Array Int Bool

nullToggleOptions :: ToggleOptions
nullToggleOptions = array (0, optionNumber - 1) (zip [0 .. optionNumber - 1] (repeat False))

filterNumber :: Int
filterNumber = 4

optionNumber :: Int
optionNumber = 8

[filterMassDensity, filterChargeDensity, filterElectronegativity, filterElectronDensity, continuousBondColoring, showLonePairs, showEmptyOrbitals, showConvexHulls] = [0 .. optionNumber - 1]

data RendAtom = AtomSphere Vector Double ColorRGBA
  deriving (Show)

instance Transform3D RendAtom where
  translate3D v     (AtomSphere p r col) = AtomSphere (translate3D v     p) r col
  rotate3D axis ang (AtomSphere p r col) = AtomSphere (rotate3D axis ang p) r col
  scale3D k         (AtomSphere p r col) = AtomSphere (scale3D k         p) (r * k) col

initLonePair :: ColorRGBA -> [RendAtom]
initLonePair (ColorRGBA r g b a) = [ AtomSphere defaultLoneElectronVector1 defaultLoneElectronRadius defaultLoneElectronColor
                                   , AtomSphere defaultLoneElectronVector2 defaultLoneElectronRadius defaultLoneElectronColor
                                   , AtomSphere origin defaultLonePairRadius (ColorRGBA r g b (a * defaultLonePairTransparency)) ]

data RendBond = BondCylinder (Maybe (Double,Double)) Vector Double ColorRGBA Vector Double ColorRGBA Double
  deriving (Show)

instance Transform3D RendBond where
  translate3D v     (BondCylinder jrr s sl sCol e el eCol r) = BondCylinder jrr (translate3D v     s) sl sCol (translate3D v     e) el eCol r
  rotate3D axis ang (BondCylinder jrr s sl sCol e el eCol r) = BondCylinder jrr (rotate3D axis ang s) sl sCol (rotate3D axis ang e) el eCol r
  scale3D k (BondCylinder jrr s sl sCol e el eCol r) = BondCylinder (fmap (\(a,b) -> (a * k, b * k)) jrr) (scale3D k s) (sl * k) sCol (scale3D k e) (el * k) eCol (r * k)


data RendHullVertex = HullSphere Point Double
  deriving (Show)

instance Transform3D RendHullVertex where
  translate3D v     (HullSphere p r) = HullSphere (translate3D v     p) r
  rotate3D axis ang (HullSphere p r) = HullSphere (rotate3D axis ang p) r
  scale3D k         (HullSphere p r) = HullSphere (scale3D k         p) (r * k)


data RendHullEdge = HullCylinder Point Double Point Double
  deriving (Show)

instance Transform3D RendHullEdge where
  translate3D v     (HullCylinder s sr e er) = HullCylinder (translate3D v     s) sr (translate3D v     e) er
  rotate3D axis ang (HullCylinder s sr e er) = HullCylinder (rotate3D axis ang s) sr (rotate3D axis ang e) er
  scale3D k (HullCylinder s sr e er) = HullCylinder (scale3D k s) (sr * k) (scale3D k e) (er * k)


data RendHullPoly = HullPolygon Polygon Vector
  deriving (Show)

instance Transform3D RendHullPoly where
  translate3D v     (HullPolygon ps n) = HullPolygon (translate3D v     ps) n
  rotate3D axis ang (HullPolygon ps n) = HullPolygon (rotate3D axis ang ps) (rotate3D axis ang n)
  scale3D k         (HullPolygon ps n) = HullPolygon (scale3D k         ps) n


data RendHull = RendHull [RendHullVertex] [RendHullEdge] [RendHullPoly] ColorRGBA
  deriving (Show)

instance Transform3D RendHull where
  translate3D v     (RendHull rvs res rps col) = RendHull (translate3D v     rvs) (translate3D v     res) (translate3D v     rps) col
  rotate3D axis ang (RendHull rvs res rps col) = RendHull (rotate3D axis ang rvs) (rotate3D axis ang res) (rotate3D axis ang rps) col
  scale3D k         (RendHull rvs res rps col) = RendHull (scale3D k         rvs) (scale3D k         res) (scale3D k         rps) col


data RendMolecule = RendMolecule Double [RendAtom] [RendBond] [RendHull]
  deriving (Show)

nullRendMol :: RendMolecule
nullRendMol = RendMolecule 1 [] [] []

instance Transform3D RendMolecule where
  translate3D v     (RendMolecule size ras rbs rhs) = RendMolecule  size      (translate3D v     ras) (translate3D v     rbs) (translate3D v     rhs)
  rotate3D axis ang (RendMolecule size ras rbs rhs) = RendMolecule  size      (rotate3D axis ang ras) (rotate3D axis ang rbs) (rotate3D axis ang rhs)
  scale3D k         (RendMolecule size ras rbs rhs) = RendMolecule (size * k) (scale3D k         ras) (scale3D k         rbs) (scale3D k         rhs)


-- Default Rendering Parameters

defaultSingleBondRadius :: Double
defaultSingleBondRadius = 10

defaultDoubleBondRadius :: Double
defaultDoubleBondRadius = 4

defaultDoubleBondSeparation :: Double
defaultDoubleBondSeparation = 6

defaultTripleBondRadius :: Double
defaultTripleBondRadius = 3

defaultTripleBondSeparation :: Double
defaultTripleBondSeparation = 7

defaultLonePairRadius :: Double
defaultLonePairRadius = 10

defaultLonePairTransparency :: Float
defaultLonePairTransparency = 0.3

defaultLoneElectronRadius :: Double
defaultLoneElectronRadius = 2

defaultLoneElectronColor :: ColorRGBA
defaultLoneElectronColor = ColorRGBA 0 0 0 1

defaultLoneElectronVector1 :: Vector
defaultLoneElectronVector1 = defaultLonePairRadius |*| Vector 0.3 0 0.5

defaultLoneElectronVector2 :: Vector
defaultLoneElectronVector2 = defaultLonePairRadius |*| Vector (-0.3) 0 0.5

defaultEmptyOrbitalRadius :: Double
defaultEmptyOrbitalRadius = defaultLonePairRadius

defaultEmptyOrbitalTransparency :: Float
defaultEmptyOrbitalTransparency = defaultLonePairTransparency

defaultDativeBondRadius :: Double
defaultDativeBondRadius = 7

defaultDativeBondMinimumSeparation :: Double
defaultDativeBondMinimumSeparation = 5

-- Data Transformation Functions

moleculeGetRends :: ViewMode -> ToggleOptions -> ElectronegativityType -> Molecule -> RendMolecule
moleculeGetRends vm tgs ent (Molecule ci _ as'' bs'' es _ hulls _) =
  let hullVerts = map hullVertices hulls
      (has,nas) = partition ((`elem` concat hullVerts) . fst) (zip [0 .. ci - 1] as'')
      showingHull = tgs ! showConvexHulls
      (as',bs') = if showingHull
                   then (map snd nas, filter (bondOutsideHull hullVerts) bs'')
                   else (as'', bs'')
      (as,bs) = case vm of
                  Skeletal -> (as', bs') -- (filter ((/= H) . element) as', filter (notHBond as') bs')
                  _        -> (as', bs')
      fs = getFilters (take filterNumber $ assocs tgs)
      filtered = not (null fs)
      cs = if filtered
            then map sum (transpose $ map (filterColor as) fs)
            else repeat undefined
      (ras,ds1)         = unzip  (map (toRA       filtered) (zip as cs))
      (xas,ds2,xbs,xes) = unzip4 (map (toRB  as'' filtered cs) bs)
      (eas,ds3,ebs)     = unzip3 (map (toEnd as'' filtered cs) (es ++ concat xes))
      (rhs,dsh) = unzip $ if showingHull
                           then map (toRendHull has) hulls
                           else []
  in RendMolecule (maximum (1 : ds1 ++ concat (ds2 ++ ds3 ++ dsh))) (ras ++ concat xas ++ concat eas) (concat xbs ++ concat ebs) rhs
  where bondOutsideHull vss (Bond (a,_) (b,_)) = all (\vs -> a `notElem` vs || b `notElem` vs) vss
        bondOutsideHull vss (IMF a b) = all (\vs -> a `notElem` vs || b `notElem` vs) vss

        notHBond as (Bond (x,_) (y,_)) = (element (as !! x) /= H) && (element (as !! y) /= H)
        notHBond as (IMF x y)          = (element (as !! x) /= H) && (element (as !! y) /= H)

        getFilters ((i,b):ibs) = if b
                                  then i : getFilters ibs
                                  else getFilters ibs
        getFilters _ = []

        filterColor as 0 =
          let xs = map (realToFrac . toMassDensity) as
          in map (toColor (minimum xs , maximum xs)) xs
            where toMassDensity (Atom {element = e}) = (molarMasses ! e) / (atomicRadii ! e)^3
                  toColor itv x = let b = transformToInterval (0.3,0.8) itv x
                                  in ColorRGBA 0 0 b 1
        filterColor as 1 =
          let xs = map (realToFrac . toChargeDensity) as
          in map (toColor (minimum xs , maximum xs)) xs
             where toChargeDensity (Atom {element = e, charge = q}) = fromIntegral q / (atomicRadii ! e)^3
                   toColor itv x = let g = transformToInterval (0.3,0.8) itv x
                                   in ColorRGBA 0 g 0 1
        filterColor as 2 =
          let xs = map (realToFrac . toElectronegativity) as
          in map (toColor (minimum xs , maximum xs)) xs
             where toElectronegativity (Atom {element = e}) = (electronegativities ! e) ! ent
                   toColor itv x = let r = transformToInterval (0.3,0.8) itv x
                                   in ColorRGBA r 0 0 1

        colorSum cs = let ColorRGBA r g b a = sum cs
                      in  ColorRGBA (transformToInterval (0,1) (0,1) r)
                                    (transformToInterval (0,1) (0,1) g)
                                    (transformToInterval (0,1) (0,1) b) a

        toRA filtered (Atom v e _ _ _, c) =
          let r = atomicRadii ! e
              r' = case vm of
                     Default      -> r * defaultScaleRadius
                     SpaceFilling -> r
                     Skeletal     -> defaultSingleBondRadius
              col = if filtered
                     then c
                     else colors ! e
          in (AtomSphere v r' col, r' + norm v)

        toRB as filtered cs (Bond sio@(sa,si) eio@(ea,ei)) =
          let Atom v1 e1 _ _ os1 = as !! sa
              Atom v2 e2 _ _ os2 = as !! ea
              BondSite ots bt ov1 rx = os1 !! si
              BondSite _   _  ov2 zz = os2 !! ei
              l1 = norm ov1
              l2 = norm ov2
              col1 = if filtered
                      then cs !! sa
                      else colors ! e1
              col2 = if filtered
                      then cs !! ea
                      else colors ! e2
              jrr = if tgs ! continuousBondColoring
                     then case vm of
                            Default      -> Just ((atomicRadii ! e1) * defaultScaleRadius, (atomicRadii ! e2) * defaultScaleRadius)
                            SpaceFilling -> Nothing
                            Skeletal     -> Just (defaultSingleBondRadius, defaultSingleBondRadius)
                     else Nothing
          in case bt of
               Bonding -> case vm of
                            Default ->
                              case length ots of
                                1 -> ([], [], [BondCylinder jrr v1 l1 col1 v2 l2 col2 defaultSingleBondRadius], [])
                                2 -> let k = ((v2 - v1) * rx) `scaleTo` (1.5 * defaultDoubleBondSeparation)
                                     in ([], [], [ BondCylinder jrr (v1 + k) l1 col1 (v2 + k) l2 col2 (1.5 * defaultDoubleBondRadius)
                                                 , BondCylinder jrr (v1 - k) l1 col1 (v2 - k) l2 col2 (1.5 * defaultDoubleBondRadius) ], [])
                                3 -> let k = ((v2 - v1) * rx) `scaleTo` (1.5 * defaultTripleBondSeparation)
                                     in ([], [], [ BondCylinder jrr (v1 + k) l1 col1 (v2 + k) l2 col2 (1.5 * defaultTripleBondRadius)
                                                 , BondCylinder jrr  v1      l1 col1  v2      l2 col2 (1.5 * defaultTripleBondRadius)
                                                 , BondCylinder jrr (v1 - k) l1 col1 (v2 - k) l2 col2 (1.5 * defaultTripleBondRadius) ], [])
                            SpaceFilling -> ([],[],[],[])
                            Skeletal ->
                              case length ots of
                                1 -> ([], [], [BondCylinder jrr v1 l1 col1 v2 l2 col2 defaultSingleBondRadius], [])
                                2 -> let k = ((v2 - v1) * rx) `scaleTo` defaultDoubleBondSeparation
                                     in ([], [], [ BondCylinder jrr (v1 + k) l1 col1 (v2 + k) l2 col2 defaultDoubleBondRadius
                                                 , BondCylinder jrr (v1 - k) l1 col1 (v2 - k) l2 col2 defaultDoubleBondRadius ], [])
                                3 -> let k = ((v2 - v1) * rx) `scaleTo` defaultTripleBondSeparation
                                     in ([], [], [ BondCylinder jrr (v1 + k) l1 col1 (v2 + k) l2 col2 defaultTripleBondRadius
                                                 , BondCylinder jrr  v1      l1 col1  v2      l2 col2 defaultTripleBondRadius
                                                 , BondCylinder jrr (v1 - k) l1 col1 (v2 - k) l2 col2 defaultTripleBondRadius ], [])
               _ -> let (sr, er, lpr', eor', dr, dms) =
                          case vm of
                            Default -> ( (atomicRadii ! e1) * defaultScaleRadius
                                       , (atomicRadii ! e2) * defaultScaleRadius
                                       , defaultLonePairRadius
                                       , defaultEmptyOrbitalRadius
                                       , defaultDativeBondRadius
                                       , defaultDativeBondMinimumSeparation )
                            SpaceFilling -> ( atomicRadii ! e1
                                            , atomicRadii ! e2
                                            , defaultLonePairRadius              / defaultScaleRadius
                                            , defaultEmptyOrbitalRadius          / defaultScaleRadius
                                            , defaultDativeBondRadius            / defaultScaleRadius
                                            , defaultDativeBondMinimumSeparation / defaultScaleRadius )
                            Skeletal -> ( defaultSingleBondRadius
                                        , defaultSingleBondRadius
                                        , 0
                                        , 0
                                        , defaultSingleBondRadius
                                        , defaultSingleBondRadius * defaultDativeBondRadius / defaultDativeBondMinimumSeparation )
                        lpr = if tgs ! showLonePairs
                               then lpr'
                               else 0
                        eor = if tgs ! showEmptyOrbitals
                               then eor'
                               else 0
                        (sp,ep) = case bt of
                                    LonePair  -> ( v1 + (ov1 `scaleTo` (sr + lpr)), v2 + (ov2 `scaleTo` (er + eor)) )
                                    Accepting -> ( v1 + (ov1 `scaleTo`  sr)       , v2 + (ov2 `scaleTo`  er)        )
                        v = ep - sp
                        l = norm v
                        l' = l + 2 * dr
                        n = floor (l' / (2 * dr + dms))
                        space = l' / (fromIntegral n)
                        ds = case [space - dr, 2 * space - dr .. l] of
                               [] -> []
                               xs -> init xs
                        cs = if tgs ! continuousBondColoring
                              then map ((mixColor col1 col2) . realToFrac . (/ l)) ds
                              else repeat col1
                        spheres = zipWith (toSphere sp (unit v) dr) cs ds
                    in (spheres, [], [], [End origin Bonding sio, End origin Bonding eio])

        toRB as filtered cs' (IMF sa ea) =
          let Atom v1 e1 _ _ os1 = as !! sa
              Atom v2 e2 _ _ os2 = as !! ea
              col1 = if filtered
                      then cs' !! sa
                      else colors ! e1
              col2 = if filtered
                      then cs' !! ea
                      else colors ! e2
              (sr, er, dr, dms) =
                case vm of
                  Default -> ( (atomicRadii ! e1) * defaultScaleRadius
                             , (atomicRadii ! e2) * defaultScaleRadius
                             , defaultDativeBondRadius
                             , defaultDativeBondMinimumSeparation )
                  SpaceFilling -> ( atomicRadii ! e1
                                  , atomicRadii ! e2
                                  , defaultDativeBondRadius            / defaultScaleRadius
                                  , defaultDativeBondMinimumSeparation / defaultScaleRadius )
                  Skeletal -> ( defaultSingleBondRadius
                              , defaultSingleBondRadius
                              , defaultSingleBondRadius
                              , defaultSingleBondRadius * defaultDativeBondRadius / defaultDativeBondMinimumSeparation )
              sp = v1 + ((v2 - v1) `scaleTo` sr)
              ep = v2 + ((v1 - v2) `scaleTo` er)
              v = ep - sp
              l = norm v
              l' = l + 2 * dr
              n = floor (l' / (2 * dr + dms))
              space = l' / (fromIntegral n)
              ds = case [space - dr, 2 * space - dr .. l] of
                     [] -> []
                     xs -> init xs
              cs = map ((mixColor col1 col2) . realToFrac . (/ l)) ds
              spheres = zipWith (toSphere sp (unit v) dr) cs ds
          in (spheres, [], [], [])

        mixColor (ColorRGBA r1 g1 b1 a1) (ColorRGBA r2 g2 b2 a2) n =
          let m = 1 - n
              r = m * r1 + n * r2
              g = m * g1 + n * g2
              b = m * b1 + n * b2
              a = m * a1 + n * a2
          in ColorRGBA r g b a

        toSphere sp u dr col d = AtomSphere (sp + (d |*| u)) dr col

        toEnd as filtered cs (End p' _ (i,o)) =
          let Atom p e _ _ os = as !! i
              BondSite ots et v rx = os !! o
              l   = norm v
              col = if filtered
                     then cs !! i
                     else colors ! e
              rbs = case et of
                      Bonding -> case vm of
                                   SpaceFilling -> []
                                   _ -> case length ots of
                                          1 -> [BondCylinder Nothing p l col p' 0 colorless defaultSingleBondRadius]
                                          2 -> let k = (v * rx) `scaleTo` defaultDoubleBondSeparation
                                               in [ BondCylinder Nothing (p + k) l col (p' + k) 0 colorless defaultDoubleBondRadius
                                                  , BondCylinder Nothing (p - k) l col (p' - k) 0 colorless defaultDoubleBondRadius ]
                      _ -> []
              (ras,ds) = case et of
                           Bonding  -> unzip (map toCap rbs)
                           LonePair -> if tgs ! showLonePairs
                                        then case vm of
                                               Skeletal -> ([],[])
                                               _ -> let (r,scl) = case vm of
                                                                    Default      -> ((atomicRadii ! e)    * defaultScaleRadius, 1)
                                                                    SpaceFilling -> ( atomicRadii ! e , 1 / defaultScaleRadius   )
                                                        lp = map ( (translate3D p)
                                                                 . (alignTo v   (Vector 0 0 1))
                                                                 . (translate3D (Vector 0 0 r))
                                                                 . (scale3D scl)) (initLonePair col)
                                                        AtomSphere pos _ _ = last lp
                                                    in (lp, [norm pos + scl * defaultLonePairRadius])
                                        else ([],[])
                           Accepting -> if tgs ! showEmptyOrbitals
                                         then case vm of
                                                Skeletal -> ([],[])
                                                _ -> let scl = case vm of
                                                                 Default      -> defaultScaleRadius
                                                                 SpaceFilling -> 1
                                                         r = (atomicRadii ! e) * scl
                                                         pos = p + v `scaleTo` r
                                                         lpr = defaultEmptyOrbitalRadius
                                                         ColorRGBA cr cg cb ca = col
                                                     in ([AtomSphere pos lpr (ColorRGBA cr cg cb (ca * defaultEmptyOrbitalTransparency))], [norm pos + lpr])
                                         else ([],[])
          in (ras, ds, rbs)
          where toCap (BondCylinder _ _ _ col p _ _ r) = (AtomSphere p r col, r + norm p)

        toRendHull has (IndexHull is col) =
          let ConvexHull vs es fs = convexHull (nub (map toHullPoint is))
              (rvs,ds) = unzip (map toRendHullVertex vs)
          in (RendHull rvs (map toRendHullEdge es) (map toRendHullPoly fs) col, ds)
          where toHullPoint i = let Just (Atom p e _ _ _) = lookup i has
                                    r = atomicRadii ! e
                                    r' = case vm of
                                           Default      -> r * defaultScaleRadius
                                           SpaceFilling -> r
                                           Skeletal     -> defaultSingleBondRadius
                                in HullPoint p r'

                toRendHullVertex (HullPoint p r) = (HullSphere p r, norm p + r)
                toRendHullEdge (HullPoint p1 r1, HullPoint p2 r2) = HullCylinder p1 r1 p2 r2
                toRendHullPoly (HullFace hps (Plane v _)) = HullPolygon (map (toTriangle v) hps) v
                  where toTriangle v (HullPoint p r) = p + (r |*| v)

defaultQuadricStyle :: QuadricStyle
defaultQuadricStyle = QuadricStyle (Just Smooth) NoTextureCoordinates Outside FillStyle

renderMolecule :: IO () -> Double -> [Vector] -> [[Double]] -> RendMolecule -> IO ()
renderMolecule io scl ns xyz (RendMolecule _ ras rbs rhs) = do
  loadIdentity
  io
  mapM_ drawAtom ras
  mapM_ drawBond rbs
  mapM_ drawHull rhs
  where sliceAdjust r = ceiling (2.6 * sqrt (r * scl))

        inViewSpace [nx,ny,nz] [[x1,x2],[y1,y2],[z1,z2]] r pos =
          let x = nx `dot` pos
              y = ny `dot` pos
              z = nz `dot` pos
          in (  x1 - r < x && x < x2 + r
             && y1 - r < y && y < y2 + r
             && z1 - r < z && z < z2 + r )

        drawAtom (AtomSphere pos r col) = preservingMatrix $ when (inViewSpace ns xyz r pos) $ do
          let n = sliceAdjust r
          translateGL pos
          colorGL col
          renderQuadric defaultQuadricStyle (Sphere (realToFrac r) n n)

        drawBond (BondCylinder jrr s sl' sCol e el' eCol r') = preservingMatrix $ when (inViewSpace ns xyz r' s || inViewSpace ns xyz r' e) $ do
          let n  = sliceAdjust r'
              v  = e - s
              k  = norm v / (sl' + el')
              sl = k * sl'
              el = k * el'
              ang = angleBetween (Vector 0 0 1) v
              p = v * (Vector 0 0 1)
              r = realToFrac r'
          translateGL s
          if p == Vector 0 0 0
           then let k = 1 - 2 * ang / pi
                in scaleGL k
           else rotateRadGL (-ang) p
          case jrr of
            Just (sr,er) -> do
              let xys = map (\t -> (cos t, sin t)) [0, 2 * pi / fromIntegral n .. 2 * pi]
                  sz  = realToFrac sr
              colorGL sCol
              renderQuadric defaultQuadricStyle (Cylinder r r sz n 1)
              translate (Vector3 0 0 sz)
              let mz  = realToFrac (sl + el - sr - er)
                  mqs = toMidQuadStrip r sCol eCol mz xys
              renderPrimitive QuadStrip mqs
              translate (Vector3 0 0 mz)
              colorGL eCol
              renderQuadric defaultQuadricStyle (Cylinder r r (realToFrac er) n 1)
            _ -> do
              colorGL sCol
              renderQuadric defaultQuadricStyle (Cylinder r r (realToFrac sl) n 1)
              translateGL (Vector 0 0 sl)
              colorGL eCol
              renderQuadric defaultQuadricStyle (Cylinder r r (realToFrac el) n 1)
          where toMidQuadStrip r sCol eCol l ((x',y'):xys) = do
                  normal (Normal3 x' y' 0)
                  let x = r * x'
                      y = r * y'
                  colorGL sCol
                  vertex (Vertex3 x y 0)
                  colorGL eCol
                  vertex (Vertex3 x y l)
                  toMidQuadStrip r sCol eCol l xys
                toMidQuadStrip _ _ _ _ _ = return ()

        drawHull (RendHull rvs res pss col) = preservingMatrix $ do
          colorGL col
          mapM_ drawHullVertex rvs
   --       colorGL red
          mapM_ drawHullEdge   res
   --       colorGL blue
          mapM_ drawHullPoly   pss
          where drawHullVertex (HullSphere pos r) = preservingMatrix $ when (inViewSpace ns xyz r pos) $ do
                  let n = sliceAdjust r
                  translateGL pos
                  renderQuadric defaultQuadricStyle (Sphere (realToFrac r) n n)

                drawHullEdge (HullCylinder p1' r1' p2' r2') = preservingMatrix $ when (inViewSpace ns xyz r1' p1' || inViewSpace ns xyz r2' p2') $ do
                  let v' = p2' - p1'
                      (p1,p2,r1,r2) = if r1' < r2'
                                       then let k = (r2' - r1') / norm v'
                                                x1 = k * r1'
                                                x2 = k * r2'
                                            in (p1' - (v' `scaleTo` x1), p2' - (v' `scaleTo` x2), sqrt (r1'^2 - x1^2), sqrt (r2'^2 - x2^2))
                                       else let k = (r1' - r2') / norm v'
                                                x1 = k * r1'
                                                x2 = k * r2'
                                            in (p1' + (v' `scaleTo` x1), p2' + (v' `scaleTo` x2), sqrt (r1'^2 - x1^2), sqrt (r2'^2 - x2^2))
                      n = sliceAdjust (max r1 r2)
                      v = p2 - p1
                      ang = angleBetween (Vector 0 0 1) v
                      p = v * (Vector 0 0 1)
                  translateGL p1
                  if p == Vector 0 0 0
                   then let k = 1 - 2 * ang / pi
                        in scaleGL k
                   else rotateRadGL (-ang) p
                  renderQuadric defaultQuadricStyle (Cylinder (realToFrac r1) (realToFrac r2) (realToFrac (norm v)) n 1)

                drawHullPoly (HullPolygon ps v) = when (any (inViewSpace ns xyz 0) ps) $ renderPrimitive Polygon $ do
                  normalGL v
                  mapM_ vertexGL ps