{-# LANGUAGE TypeSynonymInstances #-}

module Chemistry where

import Geometry
import Data.Array

-- Atomic Data

data Element = Empty
             | H                                                                             | He        
             | Li | Be                                              | B  | C  | N  | O  | F  | Ne
             | Na | Mg                                              | Al | Si | P  | S  | Cl | Ar
             | K  | Ca | Ti | V  | Cr | Mn | Fe | Co | Ni | Cu | Zn | Ga | Ge | As | Se | Br | Kr
  deriving (Show, Eq, Ord, Enum, Read)

atomicNumber :: Element -> Int
atomicNumber = fromEnum

instance Ix Element where
  range = uncurry enumFromTo
  index   (l,_) a = fromEnum a - fromEnum l
  inRange (l,u) a = fromEnum l <= i && i <= fromEnum u
    where i = fromEnum a
  rangeSize (l,u) = fromEnum u - fromEnum l + 1

-- Van Der Waals Radii, picometer
atomicRadii :: Array Element Double
atomicRadii =
  array (H,Ar)
    [ (H , 120)
    , (He, 140)
    , (Li, 182)
    , (Be, undefined)
    , (B , undefined)
    , (C , 170)
    , (N , 155)
    , (O , 152)
    , (F , 147)
    , (Ne, 154)
    , (Na, 227)
    , (Mg, 173)
    , (Al, undefined)
    , (Si, 210)
    , (P , 180)
    , (S , 180)
    , (Cl, 175)
    , (Ar, 188) ]

roundDiameter :: Double -> Double
roundDiameter k = fromIntegral (round (k * 10^2)) / 10^2

molarMasses :: Array Element Double
molarMasses =
  array (H,Ar)
    [ (H , 1.00794)
    , (He, 4.002602)
    , (Li, 6.941)
    , (Be, 9.012182)
    , (B , 10.811)
    , (C , 12.0107)
    , (N , 14.0067)
    , (O , 15.9994)
    , (F , 18.9984032)
    , (Ne, 20.1797)
    , (Na, 22.98976928)
    , (Mg, 24.3050)
    , (Al, 26.9815386)
    , (Si, 28.0855)
    , (P , 30.973762)
    , (S , 32.065)
    , (Cl, 35.453)
    , (Ar, 39.0983) ]

roundMass :: Double -> Double
roundMass k = fromIntegral (round (k * 10^5)) / 10^5

colors :: Array Element ColorRGBA
colors =
  array (H,Ar)
    [ (H , white)
    , (He, undefined)
    , (Li, undefined)
    , (Be, undefined)
    , (B , undefined)
    , (C , black)
    , (N , blue)
    , (O , red)
    , (F , undefined)
    , (Ne, undefined)
    , (Na, purple)
    , (Mg, undefined)
    , (Al, grey)
    , (Si, grey)
    , (P , orange)
    , (S , yellow)
    , (Cl, green)
    , (Ar, undefined) ]

data Orbital = O_1s
             | O_2s | O_2px | O_2py | O_2pz
             | O_3s | O_3px | O_3py | O_3pz | O_3dz2 | O_3dxz | O_3dyz | O_3dxy | O_3dx2_y2
             | Hybridized ![Orbital]
  deriving (Show, Eq, Read)

data OrbitalType = Accepting | Bonding | LonePair
  deriving (Show, Eq, Ord, Enum, Read)

data Hybridization = H_Null | H_s | H_sp | H_sp2 | H_sp3
  deriving (Show, Eq, Enum, Ord, Read)

emptyShell = H_Null
s   = H_s
sp  = H_sp
sp2 = H_sp2
sp3 = H_sp3

type Charge = Int

showCharge :: Charge -> String
showCharge q | q == 0    = "0"
             | q <  0    = show (-q) ++ "-"
             | otherwise = show q ++ "+"

type OxidationState = Charge

data BondSite = BondSite { orbitals  :: [Orbital]
                         , orbType   :: OrbitalType
                         , direction :: Vector
                         , refAxis   :: Vector }
  deriving (Show, Eq, Read)

instance Transform3D BondSite where
  translate3D _ bs = bs
  rotate3D axis ang (BondSite os ot dir rx) = BondSite os ot (rotate3D axis ang dir) (rotate3D axis ang rx)
  scale3D k         (BondSite os ot dir rx) = BondSite os ot (scale3D k         dir) (scale3D k         rx)

valenceOrbitals :: Array Element (Array OxidationState [ (Hybridization, [BondSite]) ])
valenceOrbitals =
  array (H,Ar)
    [ (H , array (0,0) [ (0, [(s, [BondSite [O_1s] Bonding (Vector 0 0 31) (Vector 0 0 1) ]) ]) ])
    , (He, undefined)
    , (Li, undefined)
    , (Be, undefined)
    , (B , undefined)
    , (C , array (0,0)
             [ (0, [ (sp3, let k = 77 / sqrt 3
                           in [ BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding (Vector   k    k    k  ) (Vector (-1) (-1)   2 )
                              , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding (Vector   k  (-k) (-k) ) (Vector (-1)   1  (-2))
                              , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding (Vector (-k)   k  (-k) ) (Vector   1  (-1) (-2))
                              , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding (Vector (-k) (-k)   k  ) (Vector   1    1    2 ) ])

                   , (sp2, [ BondSite [Hybridized [O_2s, O_2px, O_2py] , O_2pz]  Bonding (Vector 0 (-71) 0 )      (Vector 0 0 1)
                           , BondSite [Hybridized [O_2s, O_2px, O_2py]] Bonding (Vector ( sqrt 3 * 35.5) 35.5 0 ) (Vector 0 0 1)
                           , BondSite [Hybridized [O_2s, O_2px, O_2py]] Bonding (Vector (-sqrt 3 * 35.5) 35.5 0 ) (Vector 0 0 1) ])

                   , (sp , [ BondSite [Hybridized [O_2s, O_2px], O_2py]  Bonding (Vector   71  0 0 ) (Vector 0 0 1)
                           , BondSite [Hybridized [O_2s, O_2px], O_2pz]  Bonding (Vector (-71) 0 0 ) (Vector 0 1 0) ]) ]) ])
    , (N , array (0,1)
             [ (0, [ (sp3, let k = 71 / sqrt 3
                           in [ BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  (Vector   k    k    k  ) (Vector (-1) (-1)   2 )
                              , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  (Vector   k  (-k) (-k) ) (Vector (-1)   1  (-2))
                              , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  (Vector (-k)   k  (-k) ) (Vector   1  (-1) (-2))
                              , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] LonePair (Vector (-k) (-k)   k  ) (Vector   1    1    2 ) ])

                   , (sp2, [ BondSite [Hybridized [O_2s, O_2px, O_2py] , O_2pz]  Bonding (Vector 0 (-71) 0 )       (Vector 0 0 1)
                           , BondSite [Hybridized [O_2s, O_2px, O_2py]] Bonding  (Vector ( sqrt 3 * 35.5) 35.5 0 ) (Vector 0 0 1)
                           , BondSite [Hybridized [O_2s, O_2px, O_2py]] LonePair (Vector (-sqrt 3 * 35.5) 35.5 0 ) (Vector 0 0 1) ])

                   , (sp , [ BondSite [Hybridized [O_2s, O_2px], O_2py, O_2pz]  Bonding  (Vector   71  0 0 ) (Vector 0 0 1)
                           , BondSite [Hybridized [O_2s, O_2px]              ]  LonePair (Vector (-71) 0 0 ) (Vector 0 1 0) ]) ])

             , (1, [ (sp3, let k = 71 / sqrt 3
                           in [ BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding (Vector   k    k    k  ) (Vector (-1) (-1)   2 )
                              , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding (Vector   k  (-k) (-k) ) (Vector (-1)   1  (-2))
                              , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding (Vector (-k)   k  (-k) ) (Vector   1  (-1) (-2))
                              , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding (Vector (-k) (-k)   k  ) (Vector   1    1    2 ) ])
                   , (sp2, [ BondSite [Hybridized [O_2s, O_2px, O_2py] , O_2pz] Bonding (Vector 0 (-71) 0 )       (Vector 0 0 1)
                           , BondSite [Hybridized [O_2s, O_2px, O_2py]] Bonding (Vector ( sqrt 3 * 35.5) 35.5 0 ) (Vector 0 0 1)
                           , BondSite [Hybridized [O_2s, O_2px, O_2py]] Bonding (Vector (-sqrt 3 * 35.5) 35.5 0 ) (Vector 0 0 1) ]) ]) ])
    , (O , array (-1,0)
             [ (-1, [ (sp3, let k = 62.7 / sqrt 3
                            in [ BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  (Vector   k    k    k  ) (Vector (-1) (-1)   2 )
                               , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] LonePair (Vector   k  (-k) (-k) ) (Vector (-1)   1  (-2))
                               , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] LonePair (Vector (-k)   k  (-k) ) (Vector   1  (-1) (-2))
                               , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] LonePair (Vector (-k) (-k)   k  ) (Vector   1    1    2 ) ]) ])
             , ( 0, [ (sp3, let k = 62.7 / sqrt 3
                            in [ BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  (Vector   k    k    k  ) (Vector (-1) (-1)   2 )
                               , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  (Vector   k  (-k) (-k) ) (Vector (-1)   1  (-2))
                               , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] LonePair (Vector (-k)   k  (-k) ) (Vector   1  (-1) (-2))
                               , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] LonePair (Vector (-k) (-k)   k  ) (Vector   1    1    2 ) ])

                    , (sp2, [ BondSite [Hybridized [O_2s, O_2px, O_2py], O_2pz] Bonding (Vector 0 (-60.5) 0 )         (Vector 0 0 1)
                            , BondSite [Hybridized [O_2s, O_2px, O_2py]] LonePair (Vector ( sqrt 3 * 30.25) 30.25 0 ) (Vector 0 0 1)
                            , BondSite [Hybridized [O_2s, O_2px, O_2py]] LonePair (Vector (-sqrt 3 * 30.25) 30.25 0 ) (Vector 0 0 1) ]) ]) ])
    , (F , undefined)
    , (Ne, undefined)
    , (Na, array (1,1)
             [ (1, [ (emptyShell, []) ]) ])
    , (Mg, undefined)
    , (Al, array (3,3)
             [ (3, [ (emptyShell, []) ]) ])
    , (Si, array (0,0)
             [ (0, [ (sp3, let k = 110 / sqrt 3
                           in [ BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] Bonding  (Vector   k    k    k  ) (Vector (-1) (-1)   2 )
                              , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] Bonding  (Vector   k  (-k) (-k) ) (Vector (-1)   1  (-2))
                              , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] Bonding  (Vector (-k)   k  (-k) ) (Vector   1  (-1) (-2))
                              , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] Bonding  (Vector (-k) (-k)   k  ) (Vector   1    1    2 ) ]) ]) ])
    , (P , array (0,0)
             [ (0, [ (sp3, let k = 107 / sqrt 3
                           in [ BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] Bonding  (Vector   k    k    k  ) (Vector (-1) (-1)   2 )
                              , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] Bonding  (Vector   k  (-k) (-k) ) (Vector (-1)   1  (-2))
                              , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] Bonding  (Vector (-k)   k  (-k) ) (Vector   1  (-1) (-2))
                              , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] LonePair (Vector (-k) (-k)   k  ) (Vector   1    1    2 ) ]) ]) ])
    , (S , array (0,0)
             [ (0, [ (sp3, let k = 105 / sqrt 3
                           in [ BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] Bonding  (Vector   k    k    k  ) (Vector (-1) (-1)   2 )
                              , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] Bonding  (Vector   k  (-k) (-k) ) (Vector (-1)   1  (-2))
                              , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] LonePair (Vector (-k)   k  (-k) ) (Vector   1  (-1) (-2))
                              , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]] LonePair (Vector (-k) (-k)   k  ) (Vector   1    1    2 ) ]) ]) ])
    , (Cl, array (-1,-1)
             [ (-1, [ (emptyShell, []) ]) ])
--             , (0 , [ (s, [BondSite [O_1s] Bonding (Vector 0 0 31) (Vector 0 0 1) ]) ]) ])
    , (Ar, undefined) ]

type Electronegativity = Double

data ElectronegativityType = Allen | Pauling | AllredRochow | Mulliken | Jaffe | Sanderson
  deriving (Show, Eq, Ord, Enum, Read)

instance Ix ElectronegativityType where
  range = uncurry enumFromTo
  index   (l,_) a = fromEnum a - fromEnum l
  inRange (l,u) a = fromEnum l <= i && i <= fromEnum u
    where i = fromEnum a
  rangeSize (l,u) = fromEnum u - fromEnum l + 1

electronegativities :: Array Element (Array ElectronegativityType Electronegativity)
electronegativities =
  array (H,Ar)
    [ (H , array (Allen, Sanderson) [(Allen,2.300),(Pauling,2.20),(AllredRochow,2.20),(Mulliken,3.059),(Jaffe,2.25),(Sanderson,2.59)])
    , (He, array (Allen, Sanderson) [(Allen,4.160),(Pauling,0   ),(AllredRochow,5.50),(Mulliken,0    ),(Jaffe,3.49),(Sanderson,0   )])
    , (Li, array (Allen, Sanderson) [(Allen,0.912),(Pauling,0.98),(AllredRochow,0.97),(Mulliken,1.282),(Jaffe,0.97),(Sanderson,0.89)])
    , (Be, array (Allen, Sanderson) [(Allen,1.576),(Pauling,1.57),(AllredRochow,1.47),(Mulliken,1.987),(Jaffe,1.54),(Sanderson,1.81)])
    , (B , array (Allen, Sanderson) [(Allen,2.051),(Pauling,2.04),(AllredRochow,2.01),(Mulliken,1.828),(Jaffe,2.04),(Sanderson,2.28)])
    , (C , array (Allen, Sanderson) [(Allen,2.544),(Pauling,2.55),(AllredRochow,2.50),(Mulliken,2.671),(Jaffe,2.48),(Sanderson,2.75)])
    , (N , array (Allen, Sanderson) [(Allen,3.066),(Pauling,3.04),(AllredRochow,3.07),(Mulliken,3.083),(Jaffe,2.90),(Sanderson,3.19)])
    , (O , array (Allen, Sanderson) [(Allen,3.160),(Pauling,3.44),(AllredRochow,3.50),(Mulliken,3.125),(Jaffe,3.41),(Sanderson,3.65)])
    , (F , array (Allen, Sanderson) [(Allen,4.193),(Pauling,3.98),(AllredRochow,4.10),(Mulliken,4.438),(Jaffe,3.91),(Sanderson,4.00)])
    , (Ne, array (Allen, Sanderson) [(Allen,4.787),(Pauling,0   ),(AllredRochow,4.84),(Mulliken,4.597),(Jaffe,3.98),(Sanderson,4.50)])
    , (Na, array (Allen, Sanderson) [(Allen,0.869),(Pauling,0.93),(AllredRochow,1.01),(Mulliken,1.212),(Jaffe,0.91),(Sanderson,0.56)])
    , (Mg, array (Allen, Sanderson) [(Allen,1.293),(Pauling,1.31),(AllredRochow,1.23),(Mulliken,1.630),(Jaffe,1.37),(Sanderson,1.32)])
    , (Al, array (Allen, Sanderson) [(Allen,1.613),(Pauling,1.61),(AllredRochow,1.47),(Mulliken,1.373),(Jaffe,1.83),(Sanderson,1.71)])
    , (Si, array (Allen, Sanderson) [(Allen,1.916),(Pauling,1.90),(AllredRochow,1.74),(Mulliken,2.033),(Jaffe,2.28),(Sanderson,2.14)])
    , (P , array (Allen, Sanderson) [(Allen,2.253),(Pauling,2.19),(AllredRochow,2.06),(Mulliken,2.394),(Jaffe,2.30),(Sanderson,2.52)])
    , (S , array (Allen, Sanderson) [(Allen,2.589),(Pauling,2.58),(AllredRochow,2.44),(Mulliken,2.651),(Jaffe,2.69),(Sanderson,2.96)])
    , (Cl, array (Allen, Sanderson) [(Allen,2.869),(Pauling,3.16),(AllredRochow,2.83),(Mulliken,3.535),(Jaffe,3.10),(Sanderson,3.48)])
    , (Ar, array (Allen, Sanderson) [(Allen,3.242),(Pauling,0   ),(AllredRochow,3.20),(Mulliken,3.359),(Jaffe,3.19),(Sanderson,3.31)]) ]

-- Higher-Level Data Structures

data Atom = Atom { atomPos       :: Point
                 , element       :: Element
                 , charge        :: Charge
                 , hybridization :: Hybridization
                 , bondSites     :: [BondSite] }
  deriving (Show, Eq, Read)

instance Transform3D Atom where
  translate3D v     (Atom p e q h os) = Atom (translate3D v     p) e q h os
  rotate3D axis ang (Atom p e q h os) = Atom (rotate3D axis ang p) e q h (map (rotate3D axis ang) os)
  scale3D k         (Atom p e q h os) = Atom (scale3D k         p) e q h (map (scale3D k)         os)


class Indices a where
  addToIndices :: Int -> a -> a

instance (Indices a) => Indices [a] where
  addToIndices n = map (n `addToIndices`)


type OrbitalRef = (Int,Int)

instance Indices OrbitalRef where
  n `addToIndices` (a,b) = (a + n, b)


data Bond = Bond OrbitalRef OrbitalRef | IMF !Int !Int
  deriving (Show, Eq, Read)

instance Indices Bond where
  n `addToIndices` Bond (a,b) (c,d) = Bond (a + n, b) (c + n, d)
  n `addToIndices` IMF x y = IMF (x + n) (y + n)


data End = End Point OrbitalType OrbitalRef
  deriving (Show, Eq, Read)

instance Transform3D End where
  translate3D v     (End p a b) = End (translate3D v     p) a b
  rotate3D axis ang (End p a b) = End (rotate3D axis ang p) a b
  scale3D k         (End p a b) = End (scale3D k         p) a b

instance Indices End where
  n `addToIndices` End p ot (x,y) = End p ot (x + n, y)

bondingEnd (End _ Bonding _) = True
bondingEnd _ = False


data IMFTag = IMFTag Int
  deriving (Show, Eq, Read)

instance Indices IMFTag where
  n `addToIndices` IMFTag i = IMFTag (i + n)


data IndexHull = IndexHull {hullVertices :: [Int], hullColor :: ColorRGBA}
  deriving (Show, Eq, Read)

instance Indices IndexHull where
  n `addToIndices` IndexHull is col = IndexHull (map (+ n) is) col


data MolData = MolData {massData :: Double, chargeData :: Charge, specialData :: [SpecialData]}
  deriving (Show, Eq, Read)


data SpecialData = CrystalData { dimension :: Int, density :: Double, latticeParameters :: [Maybe Double]}
                 | PolymerData {polyLength :: Int, density :: Double}
  deriving (Show, Eq, Read)


data Molecule = Molecule { nextIndex   :: Int
                         , molCentroid :: Point
                         , vertices    :: [Atom]
                         , edges       :: [Bond]
                         , ends        :: [End]
                         , imfTags     :: [IMFTag]
                         , hullIndices :: [IndexHull]
                         , molData     :: MolData }
  deriving (Show, Eq, Read)

nullMolecule :: Molecule
nullMolecule = Molecule 0 origin [] [] [] [] [] (MolData 0 0 [])

instance Transform3D Molecule where
  translate3D v     (Molecule i cen as bs es is hulls dt) = Molecule i (translate3D v     cen) (translate3D v     as) bs (translate3D v     es) is hulls dt
  rotate3D axis ang (Molecule i cen as bs es is hulls dt) = Molecule i (rotate3D axis ang cen) (rotate3D axis ang as) bs (rotate3D axis ang es) is hulls dt
  scale3D k         (Molecule i cen as bs es is hulls dt) = Molecule i (scale3D k         cen) (scale3D k         as) bs (scale3D k         es) is hulls dt

instance Indices Molecule where
  n `addToIndices` (Molecule ci cen as bs es is hulls dt) =
    Molecule (ci + n) cen as (n `addToIndices` bs)
                             (n `addToIndices` es)
                             (n `addToIndices` is)
                             (n `addToIndices` hulls) dt

-- Rendering data Structures

data ColorRGBA = ColorRGBA !Float !Float !Float !Float
  deriving (Show, Eq, Read)

instance Num ColorRGBA where
  ColorRGBA r1' g1' b1' a1 + ColorRGBA r2' g2' b2' a2 =
    let [r1,g1,b1,r2,g2,b2] =
          map (\x -> if x == 1
                      then x
                      else x) [r1',g1',b1'
                              ,r2',g2',b2']
        r = (r1 + r2)
        g = (g1 + g2)
        b = (b1 + b2)
        a = (a1 + a2)
    in ColorRGBA r g b a
  ColorRGBA r1 g1 b1 a1 - ColorRGBA r2 g2 b2 a2 =
    let r = clamp (0,1) (r1 - r2)
        g = clamp (0,1) (g1 - g2)
        b = clamp (0,1) (b1 - b2)
        a = clamp (0,1) (a1 - a2)
    in ColorRGBA r g b a
  ColorRGBA r1 g1 b1 a1 * ColorRGBA r2 g2 b2 a2 =
    let r = clamp (0,1) (r1 * r2)
        g = clamp (0,1) (g1 * g2)
        b = clamp (0,1) (b1 * b2)
        a = clamp (0,1) (a1 * a2)
    in ColorRGBA r g b a
  abs = undefined
  signum = undefined
  fromInteger 0 = ColorRGBA 0 0 0 0

colorInterval :: (Float, Float)
colorInterval = (0.3, 0.8)

colorless = ColorRGBA 0 0 0 0
white  = ColorRGBA 0.8 0.8 0.8 1
black  = ColorRGBA 0.3 0.3 0.3 1
grey   = ColorRGBA 0.6 0.6 0.6 1
red    = ColorRGBA 0.8 0   0   1
green  = ColorRGBA 0   0.8 0   1
blue   = ColorRGBA 0   0   0.8 1
yellow = ColorRGBA 0.8 0.8 0   1
orange = ColorRGBA 1   0.5 0   1
purple = ColorRGBA 0.7 0.3 0.7 1