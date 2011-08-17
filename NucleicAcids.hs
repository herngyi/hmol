module NucleicAcids where

import Geometry
import Chemistry
import MoleculeBuilder
import Molecules

import Data.Array
import Control.Monad
import Control.Parallel

-- Nucleobase Specification

data Nucleobase = Adenine | Thymine | Uracil | Guanine | Cytosine
  deriving (Show, Eq, Ord, Enum, Read)

nucleobaseToAbbrev :: Nucleobase -> Char
nucleobaseToAbbrev Adenine  = 'A'
nucleobaseToAbbrev Thymine  = 'T'
nucleobaseToAbbrev Uracil   = 'U'
nucleobaseToAbbrev Guanine  = 'G'
nucleobaseToAbbrev Cytosine = 'C'

abbrevToNucleobase :: Char -> Nucleobase
abbrevToNucleobase 'A' = Adenine
abbrevToNucleobase 'T' = Thymine
abbrevToNucleobase 'U' = Uracil
abbrevToNucleobase 'G' = Guanine
abbrevToNucleobase 'C' = Cytosine

partner :: NucleicAcidType -> Nucleobase -> Nucleobase
partner DNA Adenine = Thymine
partner RNA Adenine = Uracil
partner _ Thymine   = Adenine
partner _ Uracil    = Adenine
partner _ Guanine   = Cytosine
partner _ Cytosine  = Guanine

abbrevPartner :: NucleicAcidType -> Char -> Char
abbrevPartner DNA 'A' = 'T'
abbrevPartner RNA 'A' = 'U'
abbrevPartner _   'T' = 'A'
abbrevPartner _   'U' = 'A'
abbrevPartner _   'G' = 'C'
abbrevPartner _   'C' = 'G'

parseNucleobaseSequence :: String -> [Nucleobase]
parseNucleobaseSequence = map abbrevToNucleobase

instance Ix Nucleobase where
  range = uncurry enumFromTo
  index   (l,_) a = fromEnum a - fromEnum l
  inRange (l,u) a = fromEnum l <= i && i <= fromEnum u
    where i = fromEnum a
  rangeSize (l,u) = fromEnum u - fromEnum l + 1

baseMolecules :: Array Nucleobase Molecule
baseMolecules =
  array (Adenine, Cytosine)
    [ (Adenine , bonding_adenine )
    , (Thymine , bonding_thymine )
    , (Uracil  , bonding_uracil  )
    , (Guanine , bonding_guanine )
    , (Cytosine, bonding_cytosine) ]

-- Nucleobase Construction

penta_flat_N_sp3 =
  let v1 =                                     Vector 71 0 0
      v2 = rotate3D (Vector 0 0 1) ( 7*pi/10) (Vector 71 0 0)
      v3 = rotate3D (Vector 0 0 1) (-7*pi/10) (Vector 71 0 0)
      v4 =                                     Vector 0 0 71
  in Molecule 1 origin [Atom origin N 0 sp3 [ BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  v1 (Vector 0 0 1)
                                            , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  v2 (Vector 0 0 1)
                                            , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  v3 (Vector 0 0 1)
                                            , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] LonePair v4 (Vector 1 0 0) ]]
       []
       [ End v1 Bonding  (0,0)
       , End v2 Bonding  (0,1)
       , End v3 Bonding  (0,2)
       , End v4 LonePair (0,3) ] [] [] (MolData 0 0 [])

hexa_flat_N_sp3 =
  let v1 = Vector 71 0 0
      v2 = Vector (-35.5) ( 35.5 * sqrt 3) 0
      v3 = Vector (-35.5) (-35.5 * sqrt 3) 0
      v4 = Vector 0 0 71
  in Molecule 1 origin [Atom origin N 0 sp3 [ BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  v1 (Vector 0 0 1)
                                            , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  v2 (Vector 0 0 1)
                                            , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding  v3 (Vector 0 0 1)
                                            , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] LonePair v4 (Vector 1 0 0) ]]
       []
       [ End v1 Bonding  (0,0)
       , End v2 Bonding  (0,1)
       , End v3 Bonding  (0,2)
       , End v4 LonePair (0,3) ] [] [] (MolData 0 0 [])

purine_N     = rotate3D (Vector 0 0 1) (-7*pi/30) penta_flat_N_sp3
pyramidine_N = rotate3D (Vector 0 0 1) (-7*pi/30)  hexa_flat_N_sp3

bonding_adenine = translate3D (Vector 285 132 0) $ beginWith purine_N (0,2) $ do
  insertMoleculeTwist penta_hexa_C_sp2 (0,1) pi
  newOrbital (1,0)
  insertMoleculeTwist penta_hexa_C_sp2 (0,0) pi
  newOrbital (2,1)
  insertMoleculeTwist penta_N_sp2 (0,2) 0
  newOrbital (3,0)
  insertMoleculeTwist penta_C_sp2 (0,0) pi
  newOrbital (1,2)
  insertMolecule (atom N sp2) (0,1)
  newOrbital (5,0)
  insertMolecule (atom C sp2) (0,0)
  newOrbital (6,2)
  insertMolecule (atom N sp2) (0,1)
  newOrbital (7,0)
  insertMolecule (atom C sp2) (0,0)
  matchEnds
  newOrbital (8,1)
  insertMoleculeTwist (atom N sp3) (0,1) (5*pi/6)
  newOrbital (9,0)
  insertMolecule (atom H s) (0,0)
  hydrogenateExcept [(0,0)]
  tagIMF 7
  tagIMF 10
  findConvexHull yellow


bonding_thymine = translate3D (Vector 415 238 0) $ beginWith pyramidine_N (0,2) $ do
  insertMoleculeTwist (atom C sp2) (0,2) pi
  newOrbital (1,1)
  insertMolecule hexa_flat_N_sp3 (0,1)
  newOrbital (2,0)
  insertMolecule (atom C sp2) (0,2)
  newOrbital (3,1)
  insertMolecule (atom C sp2) (0,1)
  newOrbital (4,0)
  insertMolecule (atom C sp2) (0,0)
  matchEnds
  newOrbital (4,2)
  insertMolecule (atom C sp3) (0,0)
  newOrbital (1,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (2,2)
  insertMolecule (atom H s) (0,0)
  newOrbital (3,0)
  insertMolecule (atom O sp2) (0,0)
  hydrogenateExcept [(0,0)]
  tagIMF 8
  tagIMF 9
  findConvexHull green

bonding_uracil = translate3D (Vector 415 238 0) $ beginWith pyramidine_N (0,2) $ do
  insertMoleculeTwist (atom C sp2) (0,2) pi
  newOrbital (1,1)
  insertMolecule hexa_flat_N_sp3 (0,1)
  newOrbital (2,0)
  insertMolecule (atom C sp2) (0,2)
  newOrbital (3,1)
  insertMolecule (atom C sp2) (0,1)
  newOrbital (4,0)
  insertMolecule (atom C sp2) (0,0)
  matchEnds
  newOrbital (4,2)
  insertMolecule (atom C sp3) (0,0)
  newOrbital (1,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (2,2)
  insertMolecule (atom H s) (0,0)
  newOrbital (3,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (4,2)
  insertMolecule (atom C sp3) (0,0)
  hydrogenateExcept [(0,0)]
  tagIMF 8
  tagIMF 9
  findConvexHull purple

bonding_guanine = translate3D (Vector 266 57 0) $ beginWith purine_N (0,2) $ do
  insertMoleculeTwist penta_hexa_C_sp2 (0,1) pi
  newOrbital (1,0)
  insertMoleculeTwist penta_hexa_C_sp2 (0,0) pi
  newOrbital (2,1)
  insertMoleculeTwist penta_N_sp2 (0,2) 0
  newOrbital (3,0)
  insertMoleculeTwist penta_C_sp2 (0,0) pi
  newOrbital (1,2)
  insertMolecule (atom N sp2) (0,1)
  newOrbital (5,0)
  insertMolecule (atom C sp2) (0,0)
  newOrbital (6,2)
  insertMolecule hexa_flat_N_sp3 (0,1)
  newOrbital (7,0)
  insertMolecule (atom C sp2) (0,2)
  matchEnds
  newOrbital (6,1)
  insertMoleculeTwist (atom N sp3) (0,0) (-pi/6)
  newOrbital (9,1)
  insertMolecule (atom H s) (0,0)
  newOrbital (7,2)
  insertMolecule (atom H s) (0,0)
  newOrbital (8,0)
  insertMolecule (atom O sp2) (0,0)
  hydrogenateExcept [(0,0)]
  tagIMF 10
  tagIMF 11
  tagIMF 12
  findConvexHull blue

bonding_cytosine = translate3D (Vector 375 205 0) $ beginWith pyramidine_N (0,2) $ do
  insertMoleculeTwist (atom C sp2) (0,2) pi
  newOrbital (1,1)
  insertMolecule (atom N sp2) (0,1)
  newOrbital (2,0)
  insertMolecule (atom C sp2) (0,0)
  newOrbital (3,2)
  insertMolecule (atom C sp2) (0,1)
  newOrbital (4,0)
  insertMolecule (atom C sp2) (0,0)
  matchEnds
  newOrbital (3,1)
  insertMoleculeTwist (atom N sp3) (0,0) (pi/6)
  newOrbital (1,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (6,2)
  insertMolecule (atom H s) (0,0)
  hydrogenateExcept [(0,0)]
  tagIMF 7
  tagIMF 2
  tagIMF 8
  findConvexHull red

adenine  = beginWith bonding_adenine  (0,0) (insertMolecule (atom H s) (0,0))
thymine  = beginWith bonding_thymine  (0,0) (insertMolecule (atom H s) (0,0))
uracil   = beginWith bonding_uracil   (0,0) (insertMolecule (atom H s) (0,0))
guanine  = beginWith bonding_guanine  (0,0) (insertMolecule (atom H s) (0,0))
cytosine = beginWith bonding_cytosine (0,0) (insertMolecule (atom H s) (0,0))

-- Nucleotide Construction

data NucleicAcidType = DNA | RNA
  deriving (Show, Eq, Ord, Enum, Read)

bonding_deoxyribose = beginWith penta_O_sp3 (0,0) $ do
  insertMoleculeStretch penta_C_sp3 (0,0) 154
  newOrbital (1,2)
  insertMoleculeStretchTwist penta_C_sp3 (0,0) 154 (-2*pi/3)
  newOrbital (2,2)
  insertMoleculeStretchTwist penta_C_sp3 (0,0) 154 (-2*pi/3)
  newOrbital (3,2)
  insertMoleculeStretchTwist penta_C_sp3 (0,0) 154 (-2*pi/3)
  addBond (Bond (4,2) (0,1))
  newOrbital (4,1)
  insertMoleculeTwist (atom C sp3) (0,0) (-7*pi/12 + pi/24)
  newOrbital (5,3)
  insertMoleculeTwist warped_oxygen (0,0) (pi/3)
  hydrogenateExcept [(1,1),(3,3),(6,1)]

bonding_ribose = beginWith penta_O_sp3 (0,0) $ do
  insertMoleculeStretch penta_C_sp3 (0,0) 154
  newOrbital (1,2)
  insertMoleculeTwist penta_C_sp3 (0,0) (-2*pi/3)
  newOrbital (2,2)
  insertMoleculeTwist penta_C_sp3 (0,0) (-2*pi/3)
  newOrbital (3,2)
  insertMoleculeTwist penta_C_sp3 (0,0) (-2*pi/3)
  addBond (Bond (4,2) (0,1))
  newOrbital (4,1)
  insertMoleculeTwist (atom C sp3) (0,0) (-7*pi/12 + pi/24)
  newOrbital (5,3)
  insertMoleculeTwist warped_oxygen (0,0) (pi/3)
  newOrbital (2,3)
  insertMolecule (atom O sp3) (0,0)
  hydrogenateExcept [(1,1),(3,3),(6,1)]

warped_phosphorous =
  let bs = (bondSites . head . vertices) (atomForceAngle P sp3 2 3 (pi/3))
      bs' = zipWith (\(BondSite _ _ dir rx) (os, ot) -> BondSite os ot dir rx) bs
              [ ( [Hybridized [O_3s, O_3px, O_3py, O_3pz], O_3dz2], Bonding )
              , ( [Hybridized [O_3s, O_3px, O_3py, O_3pz]]        , Bonding )
              , ( [Hybridized [O_3s, O_3px, O_3py, O_3pz]]        , Bonding )
              , ( [Hybridized [O_3s, O_3px, O_3py, O_3pz]]        , Bonding ) ]
      [v1,v2,v3,v4] = map direction bs'
  in Molecule 1 origin [Atom origin P 0 sp3 bs'] [] [ End v1 Bonding (0,0)
                                             , End v2 Bonding (0,1)
                                             , End v3 Bonding (0,2)
                                             , End v4 Bonding (0,3) ] [] [] (MolData 0 0 [])

warped_oxygen = atomForceAngle O sp3 0 1 (7*pi/8)

bonding_phosphate_ion = beginWith warped_phosphorous (0,0) $ do
  insertMolecule (atom O sp2) (0,0)
  newOrbital (0,1)
  insertMoleculeTwist (ion O (-1) sp3) (0,0) (pi/3)
  newOrbital (0,2)
  insertMoleculeTwist (atom O sp3) (0,0) (pi/12 - pi/24)

pentose_sugar :: NucleicAcidType -> Molecule
pentose_sugar DNA = bonding_deoxyribose
pentose_sugar RNA = bonding_ribose

nucleotide :: NucleicAcidType -> Nucleobase -> Molecule
nucleotide nat bs =
  let mol = beginWith (pentose_sugar nat) (6,1) $ do
        insertMoleculeTwist bonding_phosphate_ion (0,3) (5*pi/8 - pi/48)
        findConvexHull black
  in beginWithNotCentering (baseMolecules ! bs) (0,0) (insertMoleculeTwist mol (1,1) (pi -2*pi/3 - pi/144))

antinucleotide :: NucleicAcidType -> Nucleobase -> Molecule
antinucleotide nat nb = rotate3D (Vector 0 1 0) pi mol
  where mol = nucleotide nat (partner nat nb)

instance Ix NucleicAcidType where
  range = uncurry enumFromTo
  index   (l,_) a = fromEnum a - fromEnum l
  inRange (l,u) a = fromEnum l <= i && i <= fromEnum u
    where i = fromEnum a
  rangeSize (l,u) = fromEnum u - fromEnum l + 1

-- Nucleotide Polymerization

nucleotides :: Array (NucleicAcidType, Nucleobase) Molecule
nucleotides =
  array ((DNA, Adenine), (RNA, Cytosine))
    [ ((DNA, Adenine ), nucleotide DNA Adenine )
    , ((DNA, Thymine ), nucleotide DNA Thymine )
    , ((DNA, Uracil  ), nucleotide DNA Uracil  )
    , ((DNA, Guanine ), nucleotide DNA Guanine )
    , ((DNA, Cytosine), nucleotide DNA Cytosine)
    , ((RNA, Adenine ), nucleotide RNA Adenine )
    , ((RNA, Thymine ), nucleotide RNA Thymine )
    , ((RNA, Uracil  ), nucleotide RNA Uracil  )
    , ((RNA, Guanine ), nucleotide RNA Guanine )
    , ((RNA, Cytosine), nucleotide RNA Cytosine) ]

antinucleotides :: Array (NucleicAcidType, Nucleobase) Molecule
antinucleotides =
  array ((DNA, Adenine), (RNA, Cytosine))
    [ ((DNA, Adenine ), antinucleotide DNA Adenine )
    , ((DNA, Thymine ), antinucleotide DNA Thymine )
    , ((DNA, Uracil  ), antinucleotide DNA Uracil  )
    , ((DNA, Guanine ), antinucleotide DNA Guanine )
    , ((DNA, Cytosine), antinucleotide DNA Cytosine)
    , ((RNA, Adenine ), antinucleotide RNA Adenine )
    , ((RNA, Thymine ), antinucleotide RNA Thymine )
    , ((RNA, Uracil  ), antinucleotide RNA Uracil  )
    , ((RNA, Guanine ), antinucleotide RNA Guanine )
    , ((RNA, Cytosine), antinucleotide RNA Cytosine) ]

nucleotideJoinTo :: Molecule -> Molecule -> Molecule
nucleotideJoinTo nt1 nt2 =
  let [End _ _ io3',_] = filter bondingEnd (ends nt1)
      [_,End _ _ io5'] = filter bondingEnd (ends nt2)
      Molecule i cen as bs es is hulls dt = beginWithNotCentering nt1 io3' (insertMoleculeTwist nt2 io5' (-pi/3))
  in  Molecule i cen as bs (reverse es) is hulls dt

antinucleotideJoinTo :: Molecule -> Molecule -> Molecule
antinucleotideJoinTo nt1 nt2 =
  let [_,End _ _ io5'] = filter bondingEnd (ends nt1)
      [End _ _ io3',_] = filter bondingEnd (ends nt2)
  in beginWithNotCentering nt1 io5' (insertMoleculeTwist nt2 io3' (-pi/3))

sequenceToSingleHelix :: NucleicAcidType -> String -> Molecule
sequenceToSingleHelix nat str =
  let nb:nbs = map ((,) nat) (parseNucleobaseSequence str)
  in foldl nucleotideJoinTo (nucleotides ! nb) (map (nucleotides !) nbs)

sequenceToDoubleHelix :: NucleicAcidType -> String -> Molecule
sequenceToDoubleHelix nat str =
  let nb:nbs = map ((,) nat) (parseNucleobaseSequence str)
      helix1 = foldl     nucleotideJoinTo (    nucleotides ! nb) (map (    nucleotides !) nbs)
      helix2 = foldl antinucleotideJoinTo (antinucleotides ! nb) (map (antinucleotides !) nbs)
      imfs1  = imfTags helix1
      x:xs   = imfTags helix2
      imfs2  = x:xs
  in helix1 `par` helix2 `par` beginWith helix1 (0,0) $ do
                                 addMolecule helix2
                                 zipWithM_ (basePairIMF (nextIndex helix1)) imfs1 imfs2
  where basePairIMF n (IMFTag i) (IMFTag j) = addIMF i (j + n)

-- Sample Sequences

tata_box = sequenceToDoubleHelix DNA "TATAAA"

telomere n = sequenceToDoubleHelix DNA (take n (cycle "TTAGGG"))


transcribeToSingleHelix :: String -> Molecule
transcribeToSingleHelix str = sequenceToSingleHelix RNA (map (abbrevPartner RNA) str)

transcribeToDoubleHelix :: String -> Molecule
transcribeToDoubleHelix str = sequenceToDoubleHelix RNA (map (abbrevPartner RNA) str)