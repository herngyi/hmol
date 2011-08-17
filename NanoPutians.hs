module NanoPutians where

import Geometry
import Chemistry
import MoleculeBuilder
import Molecules

data NanoPutian = NanoAthlete | NanoPilgrim | NanoGreenBeret | NanoJester
                | NanoMonarch | NanoTexan   | NanoScholar    | NanoBaker
                | NanoChef    | NanoKid     | NanoBalletDancer
  deriving (Show, Eq, Ord, Enum)



nanoBalletDancer_torso = beginWith aromatic_ring (4,1) $ do
  insertMoleculeTwist hexa_C_sp3 (0,0) (-pi/2)
  newOrbital (6,3)
  insertMoleculeTwist aromatic_ring (0,1) (-pi/2)

nanoPutianLeg = beginWith kyne (1,1) $ do
  insertMolecule (alkane 3) (0,0)

nanoPutianArm = beginWith kyne (1,1) $ do
  insertMolecule (atom C sp3) (0,0)
  newOrbital (2,1)
  insertMolecule (atom C sp3) (0,0)
  newOrbital (2,2)
  insertMolecule (atom C sp3) (0,0)
  newOrbital (2,3)
  insertMolecule (atom C sp3) (0,0)

hexa_flat_C_sp3 =
  let v1 = Vector 77 0 0
      v2 = Vector (-38.5) ( 38.5 * sqrt 3) 0
      v3 = Vector (-38.5) (-38.5 * sqrt 3) 0
      v4 = Vector 0 0 77
  in Molecule 1 origin [Atom origin C 0 sp3 [ BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding v1 (Vector 0 0 1)
                                            , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding v2 (Vector 0 0 1)
                                            , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding v3 (Vector 0 0 1)
                                            , BondSite [Hybridized [O_2s, O_2px, O_2py, O_2pz]] Bonding v4 (Vector 1 0 0) ]]
       []
       [ End v1 Bonding (0,0)
       , End v2 Bonding (0,1)
       , End v3 Bonding (0,2)
       , End v4 Bonding (0,3) ] [] [] (MolData 0 0 [])

nanoBalletDancer_head = beginWith hexa_flat_C_sp3 (0,1) $ do
  insertMoleculeStretchTwist hexa_O_sp3 (0,0) 154 (-pi/6)
  newOrbital (0,2)
  insertMoleculeStretchTwist hexa_O_sp3 (0,0) 154 (5*pi/6)
  newOrbital (1,1)
  insertMoleculeTwist hexa_C_sp3 (0,0) (-2*pi/3)
  newOrbital (2,1)
  insertMoleculeTwist hexa_C_sp3 (0,0) (-2*pi/3)
  newOrbital (3,3)
  insertMolecule hexa_C_sp3 (0,0)
  matchEnds
  newOrbital (5,1)
  insertMolecule (atom C sp3) (0,0)
  newOrbital (5,2)
  insertMolecule (atom C sp3) (0,0)

nanoBalletDancer = beginWith nanoBalletDancer_torso (0,1) $ do
  insertMoleculeTwist nanoPutianLeg (0,1) (pi/2)
  newOrbital (2,1)
  insertMoleculeTwist nanoPutianLeg (0,1) (-pi/2)
  newOrbital (9,1)
  insertMoleculeTwist nanoPutianArm (0,1) (-pi/2)
  newOrbital (11,1)
  insertMoleculeTwist nanoPutianArm (0,1) (pi/2)
  newOrbital (10,1)
  insertMolecule nanoBalletDancer_head (0,0)
  hydrogenate