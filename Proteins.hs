module Proteins where

import Geometry
import Chemistry
import Molecules
import MoleculeBuilder

import Data.Array

-- Amino Acid Specification and Abbreviations

data AminoAcid = Alanine      | Arginine  | Asparagine | AsparticAcid  | Cysteine
               | GlutamicAcid | Glutamine | Glycine    | Histidine     | Isoleucine
               | Leucine      | Lysine    | Methionine | Phenylalanine | Proline
               | Serine       | Threonine | Tryptophan | Tyrosine      | Valine
  deriving (Show, Eq, Ord, Enum, Read)

data AminoAbbrev = Ala | Arg | Asn | Asp | Cys
                 | Glu | Gln | Gly | His | Ile
                 | Leu | Lys | Met | Phe | Pro
                 | Ser | Thr | Trp | Tyr | Val
  deriving (Show, Eq, Ord, Enum, Read)

abbrevToAminoAcid :: AminoAbbrev -> AminoAcid
abbrevToAminoAcid = toEnum . fromEnum

aminoAcidToAbbrev :: AminoAcid -> AminoAbbrev
aminoAcidToAbbrev = toEnum . fromEnum

-- Variable R Groups of Amino Acids

instance Ix AminoAcid where
  range = uncurry enumFromTo
  index   (l,_) a = fromEnum a - fromEnum l
  inRange (l,u) a = fromEnum l <= i && i <= fromEnum u
    where i = fromEnum a
  rangeSize (l,u) = fromEnum u - fromEnum l + 1

aminoVariableGroups :: Array AminoAcid (Molecule, Double)
aminoVariableGroups =
  array (Alanine, Valine)
    [ (Alanine       , (methyl      , pi/3 ))
    , (Arginine      , (n_propylamidineamine, pi/3 ))
    , (Asparagine    , (acetamide   , pi/3 ))
    , (AsparticAcid  , (ethanoic    , pi/3 ))
    , (Cysteine      , (methanethiol, pi/3 ))
    , (GlutamicAcid  , (propanoic   , pi/3 ))
    , (Glutamine     , (propanamide , pi/3 ))
    , (Glycine       , (atom H s    , 0    ))
    , (Histidine     , (methylimidazole , pi/3 ))
    , (Isoleucine    , (butyl       , pi/3 ))
    , (Leucine       , (isobutyl    , pi/3 ))
    , (Lysine        , (aminobutyl  , pi/3 ))
    , (Methionine    , (ethylmethylsulfide, pi/3 ))
    , (Phenylalanine , (benzyl      , pi/3 ))
    , (Proline       ,  undefined)
    , (Serine        , (methanol    , pi/3 ))
    , (Threonine     , (ethanol     , pi/3 ))
    , (Tryptophan    , (methylindole, pi/3 ))
    , (Tyrosine      , (methylphenol, pi/3 ))
    , (Valine        , (propyl      , pi/3 )) ]


methyl = beginWith (atom C sp3) (0,0) (hydrogenateExcept [(0,0)])

n_propylamidineamine = beginWith (alkane 3) (2,3) $ do
  insertMoleculeTwist (atom N sp3) (0,0) pi
  newOrbital (9,1)
  insertMoleculeTwist (atom C sp2) (0,1) (-pi/2)
  newOrbital (10,2)
  insertMolecule (atom N sp3) (0,0)
  newOrbital (10,0)
  insertMolecule (atom N sp2) (0,0)
  hydrogenateExcept [(0,0)]

acetamide = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist (atom C sp2) (0,1) (-pi/2)
  newOrbital (1,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (1,2)
  insertMoleculeTwist (atom N sp3) (0,0) (-pi/2)
  hydrogenateExcept [(0,0)]

ethanoic = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist (atom C sp2) (0,1) pi
  newOrbital (1,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (1,2)
  insertMoleculeTwist (atom O sp3) (0,0) (-pi/6)
  hydrogenateExcept [(0,0)]

methanethiol = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist (atom S sp3) (0,0) (pi/3)
  hydrogenateExcept [(0,0)]

propanoic = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (1,3)
  insertMoleculeTwist (atom C sp2) (0,1) pi
  newOrbital (2,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (2,2)
  insertMoleculeTwist (atom O sp3) (0,0) (-pi/6)
  hydrogenateExcept [(0,0)]

propanamide = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (1,3)
  insertMoleculeTwist (atom C sp2) (0,1) pi
  newOrbital (2,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (2,2)
  insertMoleculeTwist (atom N sp3) (0,0) pi
  hydrogenateExcept [(0,0)]

-- R(Glycine) = H

methylimidazole = beginWith (atom C sp3) (0,3) $ do
  insertMolecule penta_C_sp2 (0,1)
  newOrbital (1,0)
  insertMoleculeTwist penta_C_sp2 (0,0) pi
  newOrbital (2,2)
  insertMoleculeTwist (atomForceAngle N sp2 0 1 (3*pi/5)) (0,1) 0
  newOrbital (1,2)
  insertMoleculeTwist (atomForceAngle N sp3 0 1 (3*pi/5)) (0,0) (-pi/6)
  newOrbital (3,0)
  insertMolecule penta_C_sp2 (0,0)
  matchEnds
  hydrogenateExcept [(0,0)]

butyl = beginWith (atom C sp3) (0,2) $ do
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (0,3)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (2,3)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  hydrogenateExcept [(0,0)]

isobutyl = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (1,1)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (1,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  hydrogenateExcept [(0,0)]

aminobutyl = beginWith (alkane 4) (3,3) $ do
  insertMoleculeTwist (atom N sp3) (0,0) pi
  hydrogenateExcept [(0,0)]

ethylmethylsulfide = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (1,3)
  insertMoleculeTwist (atom S sp3) (0,0) (pi/3)
  newOrbital (2,1)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  hydrogenateExcept [(0,0)]

benzyl = beginWith (atom C sp3) (0,3) $ do
  insertMolecule phenyl (0,1)
  hydrogenateExcept [(0,0)]

-- Exception: Proline

methanol = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist (atom O sp3) (0,0) (pi/3)
  hydrogenateExcept [(0,0)]

ethanol = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist (atom O sp3) (0,0) 0
  newOrbital (0,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  hydrogenateExcept [(0,0)]

methylindole = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist penta_C_sp2 (0,1) pi
  newOrbital (1,2)
  insertMoleculeTwist penta_hexa_C_sp2 (0,1) 0
  newOrbital (2,0)
  insertMoleculeTwist penta_hexa_C_sp2 (0,0) pi
  newOrbital (3,1)
  insertMoleculeTwist penta_N_sp3 (0,0) (pi/6)
  newOrbital (4,2)
  insertMoleculeTwist penta_C_sp2 (0,2) (-5*pi/6)
  newOrbital (2,2)
  insertMoleculeTwist (atom C sp2) (0,2) pi
  newOrbital (6,0)
  insertMolecule (atom C sp2) (0,0)
  newOrbital (7,1)
  insertMolecule (atom C sp2) (0,2)
  newOrbital (8,0)
  insertMolecule (atom C sp2) (0,0)
  matchEnds
  hydrogenateExcept [(0,0)]

methylphenol = beginWith (atom C sp3) (0,3) $ do
  insertMolecule aromatic_ring (0,1)
  newOrbital (4,1)
  insertMoleculeTwist (atom O sp3) (0,0) (-pi/6)
  hydrogenateExcept [(0,0)]

propyl = beginWith (atom C sp3) (0,3) $ do
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (0,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  hydrogenateExcept [(0,0)]

-- Amino Acid Construction

data AminoChirality = L | D
  deriving (Show, Eq, Ord, Enum, Read)

basalAminoZwitterion :: AminoChirality -> Molecule
basalAminoZwitterion ac =
  let (a,b) = case ac of
                L -> (1,2)
                D -> (2,1)
  in beginWith (atom C sp3) (0,a) $ do
       insertMoleculeTwist (atom C sp2) (0,1) (pi/3)
       newOrbital (0,b)
       insertMoleculeTwist (ion N 1 sp3) (0,0) pi
       newOrbital (1,0)
       insertMolecule (atom O sp2) (0,0)
       newOrbital (1,2)
       insertMoleculeTwist (ion O (-1) sp3) (0,0) (-pi/6)
       hydrogenateExcept [(0,0)]

basalAminoAcid :: AminoChirality -> Molecule
basalAminoAcid ac = beginWith (bondingBasalAminoAcid ac) (1,2) $ do
  insertMoleculeTwist (atom O sp3) (0,0) (-pi/6)
  hydrogenateExcept [(0,0)]

bondingBasalAminoAcid :: AminoChirality -> Molecule
bondingBasalAminoAcid ac =
  let (a,b) = case ac of
                L -> (1,2)
                D -> (2,1)
  in beginWith (atom C sp3) (0,a) $ do
       insertMoleculeTwist (atom C sp2) (0,1) (pi/3)
       newOrbital (0,b)
       insertMoleculeTwist (atom N sp3) (0,0) pi
       newOrbital (1,0)
       insertMolecule (atom O sp2) (0,0)
       hydrogenateExcept [(0,0),(1,2),(2,1)]
       findConvexHull green

aminoAcid :: AminoChirality -> AminoAcid -> Molecule
aminoAcid L Proline = proline
aminoAcid ac am = beginWith (basalAminoAcid ac) (0,0) $ do
  let (rGroup, twist) = aminoVariableGroups ! am
      End _ _ io = head $ filter bondingEnd (ends rGroup)
  insertMoleculeTwist rGroup io twist

aminoZwitterion :: AminoChirality -> AminoAcid -> Molecule
aminoZwitterion L Proline = zwitterion_L_proline
aminoZwitterion ac am = beginWith (basalAminoZwitterion ac) (0,0) $ do
  let (rGroup, twist) = aminoVariableGroups ! am
      End _ _ io = head $ filter bondingEnd (ends rGroup)
  insertMoleculeTwist rGroup io twist

bondingAminoAcid :: AminoChirality -> AminoAcid -> Molecule
bondingAminoAcid L Proline = bonding_L_proline
bondingAminoAcid ac am = beginWith (bondingBasalAminoAcid ac) (0,0) $ do
  let (rGroup, twist) = aminoVariableGroups ! am
      End _ _ io = head $ filter bondingEnd (ends rGroup)
  insertMoleculeTwist rGroup io twist

bonding_L_proline = beginWith penta_C_sp3 (0,1) $ do
  insertMoleculeTwist (atom C sp2) (0,1) (pi/3)
  newOrbital (0,2)
  insertMoleculeStretchTwist (atomForceAngle N sp3 0 2 (3*pi/5)) (0,0) 154 (-2*pi/3)
  newOrbital (0,0)
  insertMoleculeTwist penta_C_sp3 (0,0) (-2*pi/3)
  newOrbital (3,2)
  insertMoleculeTwist penta_C_sp3 (0,0) (-2*pi/3)
  newOrbital (4,2)
  insertMoleculeTwist penta_C_sp3 (0,0) (-2*pi/3)
  addBond (Bond (5,2) (2,2))
  newOrbital (1,0)
  insertMolecule (atom O sp2) (0,0)
  hydrogenateExcept [(1,2),(2,1)]
  findConvexHullFull [0,1,2,6] green

proline = beginWith bonding_L_proline (1,2) $ do
  insertMoleculeTwist (atom O sp3) (0,0) (-pi/6)
  hydrogenate

zwitterion_L_proline = beginWith penta_C_sp3 (0,1) $ do
  insertMoleculeTwist (atom C sp2) (0,1) (pi/3)
  newOrbital (0,2)
  insertMoleculeStretchTwist (ionForceAngle N 1 sp3 0 2 (3*pi/5)) (0,0) 154 (-2*pi/3)
  newOrbital (0,0)
  insertMoleculeTwist penta_C_sp3 (0,0) (-2*pi/3)
  newOrbital (3,2)
  insertMoleculeTwist penta_C_sp3 (0,0) (-2*pi/3)
  newOrbital (4,2)
  insertMoleculeTwist penta_C_sp3 (0,0) (-2*pi/3)
  matchEnds
  newOrbital (1,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (1,2)
  insertMoleculeTwist (ion O (-1) sp3) (0,0) (-pi/6)
  hydrogenate


alanine       = aminoAcid L Alanine
arginine      = aminoAcid L Arginine
asparagine    = aminoAcid L Asparagine
aspartic_acid = aminoAcid L AsparticAcid
cysteine      = aminoAcid L Cysteine
glutamic_acid = aminoAcid L GlutamicAcid
glutamine     = aminoAcid L Glutamine
glycine       = aminoAcid L Glycine
histidine     = aminoAcid L Histidine
isoleucine    = aminoAcid L Isoleucine
leucine       = aminoAcid L Leucine
lysine        = aminoAcid L Lysine
methionine    = aminoAcid L Methionine
phenylalanine = aminoAcid L Phenylalanine

serine        = aminoAcid L Serine
threonine     = aminoAcid L Threonine
tryptophan    = aminoAcid L Tryptophan
tyrosine      = aminoAcid L Tyrosine
valine        = aminoAcid L Valine

histidine_zwitterion = aminoZwitterion L Histidine

-- Amino Acid Polymerization

aminoJoinTo :: Molecule -> Molecule -> Molecule
am1 `aminoJoinTo` am2 =
  let [_, End _ _ cio] = filter bondingEnd (ends am1)
      [End _ _ nio, _] = filter bondingEnd (ends am2)
  in beginWith am1 cio (insertMoleculeTwist am2 nio (-pi/6))

proteinAlphaHelix :: [AminoAbbrev] -> Molecule
proteinAlphaHelix abvs = let ams = map ((bondingAminoAcid L) . abbrevToAminoAcid) abvs
                         in foldr aminoJoinTo (last ams) (init ams)

alphabet :: Molecule
alphabet = proteinAlphaHelix (enumFromTo Ala Val)

prot = proteinAlphaHelix (enumFromTo Ala Trp)