{-# LANGUAGE FlexibleInstances #-}

module MoleculeBuilder where

import Geometry
import Chemistry
import Algorithms

import Data.List
import Data.Array
import Data.Maybe

import Control.Monad

-- Mimic of the State Monad

type MoleculeBuilder a = Builder Molecule OrbitalRef a

data Builder m r a = Builder {runBuilder :: m -> r -> (a, m, r)}


instance Monad (Builder Molecule OrbitalRef) where
  return a = Builder (\mol ref -> (a, mol, ref))
  f1 >>= f2 = Builder $ \mol1 ref1 ->
    let (a, mol2, ref2) = runBuilder f1 mol1 ref1
    in runBuilder (f2 a) mol2 ref2

-- Atom Specification

ionFull :: Element -> OxidationState -> Hybridization -> Maybe [(Maybe [Orbital], OrbitalType)] -> Maybe [(Int, Vector)] -> Molecule
ionFull e q h jots jfxs = Molecule 1 origin [at] [] es [] [] (MolData m q [])
  where at = case  jfxs of
               Just fxs -> warpAtom (Atom origin e q h os') fxs
               _ -> Atom origin e q h os'
        Just os = lookup h ((valenceOrbitals ! e) ! q)
        os' = case jots of
                Just ots -> zipWith (\(BondSite xs _ dir rx) (jxs,ot) -> let xs' = case jxs of
                                                                                     Just xs0 -> xs0
                                                                                     _ -> xs
                                                                         in BondSite xs' ot dir rx) os ots
                _ -> os
        es = zipWith3 End (map direction (bondSites at)) (map orbType os') (zip (repeat 0) [0 .. length os' - 1])
        m = molarMasses ! e

ionConfigWarp :: Element -> OxidationState -> Hybridization -> [OrbitalType] -> [(Int, Vector)] -> Molecule
ionConfigWarp e q h ots fxs = ionFull e q h (Just (zip (repeat Nothing) ots)) (Just fxs)

atomConfigWarp :: Element -> Hybridization -> [OrbitalType] -> [(Int, Vector)] -> Molecule
atomConfigWarp e h ots fxs = ionFull e 0 h (Just (zip (repeat Nothing) ots)) (Just fxs)

ionConfig :: Element -> OxidationState -> Hybridization -> [OrbitalType] -> Molecule
ionConfig e q h ots = ionFull e q h (Just (zip (repeat Nothing) ots)) Nothing

atomConfig :: Element -> Hybridization -> [OrbitalType] -> Molecule
atomConfig e h ots = ionFull e 0 h (Just (zip (repeat Nothing) ots)) Nothing

ionConfigFull :: Element -> OxidationState -> Hybridization -> [([Orbital], OrbitalType)] -> Molecule
ionConfigFull e q h ots = ionFull e q h (Just (map (\(xs,ot) -> (Just xs, ot)) ots)) Nothing

atomConfigFull :: Element -> Hybridization -> [([Orbital], OrbitalType)] -> Molecule
atomConfigFull e h ots = ionFull e 0 h (Just (map (\(xs,ot) -> (Just xs, ot)) ots)) Nothing

ionWarp :: Element -> OxidationState -> Hybridization -> [(Int, Vector)] -> Molecule
ionWarp e q h fxs = ionFull e q h Nothing (Just fxs)

atomWarp :: Element -> Hybridization -> [(Int, Vector)] -> Molecule
atomWarp e h fxs = ionFull e 0 h Nothing (Just fxs)

ion :: Element -> OxidationState -> Hybridization -> Molecule
ion e q h = ionFull e q h Nothing Nothing

atom :: Element -> Hybridization -> Molecule
atom e h = ionFull e 0 h Nothing Nothing

ionForceAngle :: Element -> OxidationState -> Hybridization -> Int -> Int -> Angle -> Molecule
ionForceAngle e q h a b ang =
  let av = direction (os !! a)
      bv = direction (os !! b)
      axis = unit (av * bv)
      t = (ang - angleBetween av bv) / 2
      fxs = [(a, rotate3D axis (-t) av), (b, rotate3D axis t bv)]
      at = warpAtom (Atom origin e q h os) fxs
      es = zipWith3 End (map direction (bondSites at)) (map orbType os) (zip (repeat 0) [0 .. length os - 1])
  in Molecule 1 origin [at] [] es [] [] (MolData m q [])
  where Just os = lookup h ((valenceOrbitals ! e) ! q)
        m = molarMasses ! e

atomForceAngle :: Element -> Hybridization -> Int -> Int -> Angle -> Molecule
atomForceAngle e h a b ang = ionForceAngle e 0 h a b ang

-- Monadic Molecule-Building Commands

beginWith :: Molecule -> OrbitalRef -> MoleculeBuilder () -> Molecule
beginWith mol ref mb = centerMol mol'
  where (_, mol',_) = runBuilder mb mol ref

centerMol :: Molecule -> Molecule
centerMol mol = let cen = molCentroid mol
                in translate3D (-cen) mol

beginWithNotCentering :: Molecule -> OrbitalRef -> MoleculeBuilder () -> Molecule
beginWithNotCentering mol ref mb = mol'
  where (_, mol',_) = runBuilder mb mol ref

currMolecule :: MoleculeBuilder Molecule
currMolecule = Builder (\mol ref -> (mol, mol, ref))

currOrbital :: MoleculeBuilder OrbitalRef
currOrbital = Builder (\mol ref -> (ref, mol, ref))

newMolecule :: Molecule -> MoleculeBuilder ()
newMolecule mol = Builder (\_ ref -> ((), mol, ref))

newOrbital :: OrbitalRef -> MoleculeBuilder ()
newOrbital ref = Builder (\mol _ -> ((), mol, ref))

addMolecule :: Molecule -> MoleculeBuilder ()
addMolecule mol = Builder $ \(Molecule i2 cen2 as2 bs2 es2 is2 hulls2 (MolData m2 q2 _)) ref ->
  let Molecule i1 cen1 as1 bs1 es1 is1 hulls1 (MolData m1 q1 _) = i2 `addToIndices` mol
      n2 = fromIntegral i2
      n1 = fromIntegral i1
      cen' = (1 / n1) |*| (n2 |*| cen2 + (n1 - n2) |*| cen1)
  in ((), Molecule i1 cen' (as2 ++ as1) (bs2 ++ bs1) (es2 ++ es1) (is2 ++ is1) (hulls2 ++ hulls1) (MolData (m2 + m1) (q2 + q1) []), ref)

addBond :: Bond -> MoleculeBuilder ()
addBond b@(Bond ref1 ref2) = Builder $ \(Molecule i cen as bs es is hulls dt) ref ->
  ((), Molecule i cen as (bs ++ [b]) ([ref1,ref2] `removeEndsFrom` es) is hulls dt, ref)
addBond imf = Builder $ \(Molecule i cen as bs es is hulls dt) ref -> ((), Molecule i cen as (bs ++ [imf]) es is hulls dt, ref)

removeEndsFrom :: [OrbitalRef] -> [End] -> [End]
ios `removeEndsFrom` es = filter (\(End _ _ io) -> io `notElem` ios) es

deleteBond :: Bond -> MoleculeBuilder ()
deleteBond (Bond ref1 ref2) = Builder $ \(Molecule i cen as bs es is hulls dt) ref ->
  let es' = map (newEnd as) [ref1,ref2]
  in ((), Molecule i cen as (delete (Bond ref1 ref2) bs) (es ++ es') is hulls dt, ref)
  where newEnd as io@(i,o) = let Atom p _ _ _ os = as !! i
                                 BondSite _ ot dir _ = os !! o
                             in End (p + dir) ot io
deleteBond imf = Builder $ \(Molecule i cen as bs es is hulls dt) ref -> ((), Molecule i cen as (delete imf bs) es is hulls dt, ref)

tagIMF :: Int -> MoleculeBuilder ()
tagIMF n = Builder $ \(Molecule i cen as bs es is hulls dt) ref -> ((), Molecule i cen as bs es (is ++ [IMFTag n]) hulls dt, ref)

addIMF :: Int -> Int -> MoleculeBuilder ()
addIMF m n = Builder $ \(Molecule i cen as bs es is hulls dt) ref -> ((), Molecule i cen as (bs ++ [IMF m n]) es is hulls dt, ref)

findConvexHullFull :: [Int] -> ColorRGBA -> MoleculeBuilder ()
findConvexHullFull xs col = Builder $ \(Molecule i cen as bs es is _ dt) ref -> ((), Molecule i cen as bs es is [IndexHull xs col] dt, ref)

findConvexHull :: ColorRGBA -> MoleculeBuilder ()
findConvexHull col = Builder $ \(Molecule i cen as bs es is _ dt) ref -> ((), Molecule i cen as bs es is [IndexHull [0 .. i - 1] col] dt, ref)

center :: MoleculeBuilder ()
center = Builder $ \mol ref -> ((), centerMol mol, ref)

insertMoleculeFull ::  Molecule -> OrbitalRef -> Maybe Double -> Angle -> MoleculeBuilder ()
insertMoleculeFull mol1 ref1 jl ang = Builder (\mol2 ref2 -> ((), insertMoleculeF mol1 ref1 ang mol2 ref2, ref1))
  where insertMoleculeF :: Molecule -> OrbitalRef -> Double -> Molecule -> OrbitalRef -> Molecule
        insertMoleculeF mol'@(Molecule xi1 _    as' _ _ _ _ _) (i1,o1) twist
                        mol2@(Molecule ci2 cen2 as2 bs2 es2 is2 hulls2 (MolData m2 q2 _)) io2@(i2,o2) =
          let Atom _  _  _ _ os' = as' !! i1
              Atom p2 e2 _ _ os2 = as2 !! i2
              BondSite _    _ v1 rx  = os' !! o1
              BondSite ots2 _ v2 rx2 = os2 !! o2
              axis  = unit (-v2)
              mol'' = alignTo axis v1 mol'
              rx1 = refAxis (bondSites (vertices mol'' !! i1) !! o1)
              ang = angleBetween rx2 rx1
              mol''' | ang ~= 0  = rotate3D axis twist mol''
                     | ang ~= pi = rotate3D axis twist (rotate3D axis pi mol'')
                     | otherwise = rotate3D axis twist (rotate3D (unit (rx1 * rx2)) ang mol'')
              Atom p1 e1 _ _ os1 = vertices mol''' !! i1
              trans = case jl of
                        Just l -> p2 + (v2 `scaleTo` l) - p1
                        _      -> p2 + v2 + (v2 `scaleTo` norm v1) - p1
              Molecule ci1 cen1 as1 bs1 es1 is1 hulls1 (MolData m1 q1 _) = ci2 `addToIndices` translate3D trans mol'''
              io1' = (i1 + ci2, o1)
              bs1' = Bond io2 io1' : bs1
              n2   = fromIntegral ci2
              n1   = fromIntegral xi1
              cen' = (1 / (n2 + n1)) |*| (n2 |*| cen2 + n1 |*| cen1)
          in Molecule ci1 cen' (as2 ++ as1) (bs2 ++ bs1') ([io2,io1'] `removeEndsFrom` (es2 ++ es1)) (is2 ++ is1) (hulls2 ++ hulls1) (MolData (m2 + m1) (q2 + q1) [])

insertMoleculeStretchTwist :: Molecule -> OrbitalRef -> Double -> Angle -> MoleculeBuilder ()
insertMoleculeStretchTwist mol ref l ang = insertMoleculeFull mol ref (Just l) ang

insertMolecule :: Molecule -> OrbitalRef -> MoleculeBuilder ()
insertMolecule mol ref = insertMoleculeFull mol ref Nothing 0

insertMoleculeTwist :: Molecule -> OrbitalRef -> Angle -> MoleculeBuilder ()
insertMoleculeTwist mol ref ang = insertMoleculeFull mol ref Nothing ang

insertMoleculeStretch :: Molecule -> OrbitalRef -> Double -> MoleculeBuilder ()
insertMoleculeStretch mol ref l = insertMoleculeFull mol ref (Just l) 0

-- "Smart" Commands

matchEnds :: MoleculeBuilder ()
matchEnds = Builder $ \(Molecule i cen as bs es is hulls dt) ref ->
  let (ioss, ebs) = unzip (getMatchingEnds es)
  in ((), Molecule i cen as (bs ++ ebs) (concat ioss `removeEndsFrom` es) is hulls dt, ref)
  where getMatchingEnds (e:es) = mapMaybe (getEndBond e) es ++ getMatchingEnds es
        getMatchingEnds _ = []
        getEndBond (End p1 _ io1) (End p2 _ io2) = if p1 == p2
                                                    then Just ([io1,io2], Bond io1 io2)
                                                    else Nothing

hydrogenateExcept :: [OrbitalRef] -> MoleculeBuilder ()
hydrogenateExcept ors = Builder $ \mol@(Molecule _ _ as _ es _ _ _) ref ->
  let mol' = beginWith mol (0,0) (mapM_ (addH as) (ors `removeEndsFrom` es))
  in ((), mol', ref)
  where addH as (End _ _ io@(i,o)) = do
          let BondSite ots bt _ _ = bondSites (as !! i) !! o
          when (bt == Bonding && length ots == 1) $ do
            newOrbital io
            insertMolecule (atom H s) (0,0)

hydrogenate :: MoleculeBuilder ()
hydrogenate = hydrogenateExcept []

oxygenateExcept :: [OrbitalRef] -> MoleculeBuilder ()
oxygenateExcept ors = Builder $ \mol@(Molecule _ _ as _ es _ _ _) ref ->
  let mol' = beginWith mol (0,0) (mapM_ (addO as) (ors `removeEndsFrom` es))
  in ((), mol', ref)
  where addO as (End _ _ io@(i,o)) = do
          let BondSite ots bt _ _ = bondSites (as !! i) !! o
          when (bt == Bonding && length ots == 2) $ do
            newOrbital io
            insertMolecule (atom O sp2) (0,0)

oxygenate :: MoleculeBuilder ()
oxygenate = oxygenateExcept []

-- Crystallization

convertDensity dim dsy = (dsy / 6.023) * 10**(10*dim - 23)

deleteCrystalData :: [SpecialData] -> [SpecialData]
deleteCrystalData (CrystalData _ _ _ : xs) = xs
deleteCrystalData (_:xs) = deleteCrystalData xs
deleteCrystalData _ = []

crystallize3D :: (Vector, Int) -> (Vector, Int) -> (Vector, Int) -> [( (Int,Int,Int), [Bond] )] -> MoleculeBuilder ()
crystallize3D (v1,i1) (v2,i2) (v3,i3) pbs = Builder $ \mol ref ->
  let Molecule i cen as bs es is hulls (MolData m0 q sds) = beginWith mol (0,0) $ do
        let i3s = [0 .. i3 - 1]
            i2s = [0 .. i2 - 1]
            i1s = [0 .. i1 - 1]
        mapM_ addMolecule (tail [translate3D ( (fromIntegral a |*| v1)
                                             + (fromIntegral b |*| v2)
                                             + (fromIntegral c |*| v3) ) mol | c <- i3s, b <- i2s, a <- i1s ])
        let n1 = nextIndex mol
            n2 = n1 * i1
            n3 = n2 * i2
            bonds = concat (map (crystalBonds3D n1 i1 i1s
                                                n2 i2 i2s
                                                n3 i3 i3s ) pbs)
        mapM_ addBond bonds
      m = massData (molData mol)
      sd = CrystalData 3 (convertDensity 3 $ m / abs (v1 `dot` (v2 * v3)))
                 [ Just (norm v1)
                 , Just (norm v2)
                 , Just (norm v3)
                 , Just (angleBetween v2 v3)
                 , Just (angleBetween v3 v1)
                 , Just (angleBetween v1 v2) ]
  in ((), Molecule i cen as bs es is hulls (MolData m0 q (sd: deleteCrystalData sds)) , ref)
  where crystalBonds3D n1 i1 i1s
                       n2 i2 i2s
                       n3 i3 i3s ((x,y,z), bs) = 
          concat $ mapMaybe validBond3D [let x' = x + a
                                             y' = y + b
                                             z' = z + c
                                         in ((x',y',z'), map (assign3D a b c x' y' z') bs) | c <- i3s, b <- i2s, a <- i1s ]
          where validBond3D ((x,y,z),bond) = if    0 <= x && x < i1
                                                && 0 <= y && y < i2
                                                && 0 <= z && z < i3
                                              then Just bond
                                              else Nothing
                assign3D a b c x' y' z' (Bond ref1 ref2) = Bond (locate3D a b c `addToIndices` ref1) (locate3D x' y' z' `addToIndices` ref2)
                assign3D a b c x' y' z' (IMF m n) = IMF (locate3D a b c + m) (locate3D x' y' z' + n)

                locate3D a b c = a * n1 + b * n2 + c * n3

crystallize2D :: (Vector, Int) -> (Vector, Int) -> [( (Int,Int), [Bond] )] -> MoleculeBuilder ()
crystallize2D (v1,i1) (v2,i2) pbs = Builder $ \mol ref ->
  let Molecule i cen as bs es is hulls (MolData m0 q sds) = beginWith mol (0,0) $ do
        let i2s = [0 .. i2 - 1]
            i1s = [0 .. i1 - 1]
        mapM_ addMolecule (tail [translate3D ( (fromIntegral a |*| v1)
                                             + (fromIntegral b |*| v2) ) mol | b <- i2s, a <- i1s ])
        let n1 = nextIndex mol
            n2 = n1 * i1
            bonds = concat (map (crystalBonds2D n1 i1 i1s
                                                n2 i2 i2s ) pbs)
        mapM_ addBond bonds
      m = massData (molData mol)
      sd = CrystalData 2 (convertDensity 2 $ m / norm (v1 * v2))
                 [ Just (norm v1)
                 , Just (norm v2)
                 , Nothing
                 , Nothing
                 , Nothing
                 , Just (angleBetween v1 v2) ]
  in ((), Molecule i cen as bs es is hulls (MolData m0 q (sd: deleteCrystalData sds)) , ref)
  where crystalBonds2D n1 i1 i1s
                       n2 i2 i2s ((x,y), bs) = 
          concat $ mapMaybe validBond2D [let x' = x + a
                                             y' = y + b
                                         in ((x',y'), map (assign2D a b x' y') bs) | b <- i2s, a <- i1s ]
          where validBond2D ((x,y),bonds) = if    0 <= x && x < i1
                                               && 0 <= y && y < i2
                                             then Just bonds
                                             else Nothing
                assign2D a b x' y' (Bond ref1 ref2) = Bond (locate2D a b `addToIndices` ref1) (locate2D x' y' `addToIndices` ref2)
                assign2D a b x' y' (IMF m n) = IMF (locate2D a b + m) (locate2D x' y' + n)

                locate2D a b = a * n1 + b * n2

-- Polymerization

polymerize :: Bond -> Angle -> [Bond] -> Int -> MoleculeBuilder ()
polymerize bond twist bs m = Builder $ \mol ref ->
  let mol' = beginWith mol (0,0) (addUnits (nextIndex mol) bond bs m mol)
  in ((), mol', ref)
  where addUnits _ _ _ 1 _ = return ()
        addUnits n (Bond ref1 ref2) bs m mol = do
          let ref1' = n `addToIndices` ref1
          newOrbital ref1
          insertMoleculeTwist mol ref2 twist
          bs' <- mapM thisLinkBonds bs
          addUnits n (Bond ref1' ref2) bs' (m - 1) mol
          where thisLinkBonds (Bond ref1 ref2) = do
                  let ref1' = n `addToIndices` ref1
                      ref2' = n `addToIndices` ref2
                  addBond (Bond ref1  ref2')
                  return  (Bond ref1' ref2')
                thisLinkBonds (IMF a b) = do
                  let a' = a + n
                      b' = b + n
                  addBond (IMF a  b')
                  return  (IMF a' b')

-- Molecule Projection

{-
type Coordinates = (Double, Double)

data ProjectedAtom = ProjectedAtom Int Coordinates Double
  deriving (Show)

data ProjectedBond = ProjectedBond Int Coordinates Coordinates Double
  deriving (Show)

--data ProjectedEnd = ProjectedEnd Int 
-}
