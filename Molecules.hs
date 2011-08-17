module Molecules where

import Geometry
import Chemistry
import MoleculeBuilder

import Control.Monad

hydrogen = beginWith (atom H s) (0,0) (insertMolecule (atom H s) (0,0))

water = beginWith (atom O sp3) (0,0) hydrogenate

carbon_dioxide = beginWith (atom C sp) (0,0) $ do
  insertMolecule (atom O sp2) (0,0)
  newOrbital (0,1)
  insertMolecule (atom O sp2) (0,0)

ammonia = beginWith (atom N sp3) (0,0) hydrogenate

hydrogen_peroxide = beginWith (atom O sp3) (0,0) $ do
  insertMoleculeTwist (atom O sp3) (0,0) (-pi/3)
  hydrogenate

triple_bond_C = atomConfigFull C sp [ ([Hybridized [O_2s,O_2pz], O_2px, O_2py], Bonding)
                                    , ([Hybridized [O_2s,O_2pz]              ], Bonding) ]
kyne = beginWith triple_bond_C (0,0) $ do
  insertMolecule triple_bond_C (0,0)

hydrogen_cyanide = beginWith triple_bond_C (0,0) $ do
  insertMolecule (atom N sp) (0,0)
  hydrogenate

ethene = beginWith (atom C sp2) (0,0) $ do
  insertMolecule (atom C sp2) (0,0)
  hydrogenate

alkane :: Int -> Molecule
alkane (n + 1) = beginWith (atom C sp3) (0,3) $ do
  mapM_ nextCarbon [1 .. n]
  hydrogenateExcept [(0,0),(n,3)]
  where nextCarbon m = do
          insertMoleculeTwist (atom C sp3) (0,0) pi
          newOrbital (m,3)

ethylenediamine = beginWith (atom C sp3) (0,0) $ do
  insertMolecule (atom C sp3) (0,0)
  newOrbital (0,3)
  insertMolecule (atom N sp3) (0,0)
  newOrbital (1,3)
  insertMolecule (atom N sp3) (0,0)
  hydrogenate
  findConvexHull green

nitric_acid = beginWith (ion N 1 sp2) (0,0) $ do
  insertMolecule (atom O sp2) (0,0)
  newOrbital (0,1)
  insertMoleculeTwist (atom O sp3) (0,0) (-pi/6)
  newOrbital (0,2)
  insertMoleculeTwist (atom O sp3) (0,0) (-pi/6)
  hydrogenate

sulfuric_acid = beginWith (atomConfigFull S sp3 [ ( [Hybridized [O_3s, O_3px, O_3py, O_3pz], O_3dz2]   , Bonding )
                                                , ( [Hybridized [O_3s, O_3px, O_3py, O_3pz], O_3dx2_y2], Bonding )
                                                , ( [Hybridized [O_3s, O_3px, O_3py, O_3pz]], Bonding )
                                                , ( [Hybridized [O_3s, O_3px, O_3py, O_3pz]], Bonding ) ]) (0,0) $ do

  insertMolecule (atom O sp2) (0,0)
  newOrbital (0,1)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (0,2)
  insertMoleculeTwist (atom O sp3) (0,0) (pi/3)
  newOrbital (0,3)
  insertMoleculeTwist (atom O sp3) (0,0) (pi/3)
  hydrogenate

excited_phosphorus = atomConfigFull P sp3 [ ( [Hybridized [O_3s, O_3px, O_3py, O_3pz], O_3dz2], Bonding )
                                          , ( [Hybridized [O_3s, O_3px, O_3py, O_3pz]]        , Bonding )
                                          , ( [Hybridized [O_3s, O_3px, O_3py, O_3pz]]        , Bonding )
                                          , ( [Hybridized [O_3s, O_3px, O_3py, O_3pz]]        , Bonding ) ]

phosphoric_acid = beginWith excited_phosphorus (0,0) $ do
  insertMolecule (atom O sp2) (0,0)
  newOrbital (0,1)
  insertMoleculeTwist (atom O sp3) (0,0) (pi/3)
  newOrbital (0,2)
  insertMoleculeTwist (atom O sp3) (0,0) (pi/3)
  newOrbital (0,3)
  insertMoleculeTwist (atom O sp3) (0,0) (pi/3)
  hydrogenate

tetrahedral_phosphorus = atomWarp P sp3 $ let k = 107 / sqrt 2
                                          in [ (0, Vector (-k) (-k) 0 )
                                             , (1, Vector 0 (-k) (-k) )
                                             , (2, Vector (-k) 0 (-k) ) ]

_P4 = beginWith tetrahedral_phosphorus (0,0) $ do
  insertMolecule tetrahedral_phosphorus (0,0)
  newOrbital (0,1)
  insertMoleculeTwist tetrahedral_phosphorus (0,1) (2*pi/3)
  newOrbital (0,2)
  insertMoleculeTwist tetrahedral_phosphorus (0,2) (-2*pi/3)
  matchEnds

oxygenated_tetrahedral_phosphorus =
  let k  = 107 / sqrt 2
      l  = 107 / sqrt 3
      v1 = Vector (-k) (-k) 0
      v2 = Vector 0 (-k) (-k)
      v3 = Vector (-k) 0 (-k)
      v4 = Vector   l  l   l
  in Molecule 1 origin [Atom origin P 0 sp3 [ BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]]         Bonding v1 (Vector 0 0 1)
                                            , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]]         Bonding v2 (Vector 1 0 0)
                                            , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz]]         Bonding v3 (Vector 0 1 0)
                                            , BondSite [Hybridized [O_3s, O_3px, O_3py, O_3pz], O_3dz2] Bonding v4 (Vector (-1) (-1) 2) ]]
       []
       [ End v1 Bonding (0,0)
       , End v2 Bonding (0,1)
       , End v3 Bonding (0,2)
       , End v4 Bonding (0,3) ] [] [] (MolData 0 0 [])

_P4O4 = beginWith oxygenated_tetrahedral_phosphorus (0,0) $ do
  insertMolecule  oxygenated_tetrahedral_phosphorus (0,0)
  newOrbital (0,1)
  insertMoleculeTwist oxygenated_tetrahedral_phosphorus (0,0) 0
  newOrbital (0,2)
  insertMoleculeTwist oxygenated_tetrahedral_phosphorus (0,0) 0
  matchEnds
  oxygenate

_P4O10 = beginWith excited_phosphorus (0,1) $ do
  insertMoleculeTwist (atom O sp3) (0,0) (-pi/3)
  newOrbital (0,2)
  insertMoleculeTwist (atom O sp3) (0,0) pi
  newOrbital (0,3)
  insertMoleculeTwist (atom O sp3) (0,0) (pi/3)
  newOrbital (1,1)
  insertMoleculeTwist excited_phosphorus (0,1) (-pi/3)
  newOrbital (2,1)
  insertMoleculeTwist excited_phosphorus (0,1) (-pi/3)
  newOrbital (3,1)
  insertMoleculeTwist excited_phosphorus (0,1) (-pi/3)
  newOrbital (4,2)
  insertMoleculeTwist (atom O sp3) (0,0) pi
  newOrbital (5,2)
  insertMoleculeTwist (atom O sp3) (0,0) pi
  newOrbital (6,2)
  insertMoleculeTwist (atom O sp3) (0,0) pi
  matchEnds
  oxygenate

benzene = beginWith aromatic_ring (0,1) $ do
  hydrogenate
  findConvexHullFull [0,1,2,3,4,5,6,7,8,9] green

benzoic_acid = beginWith aromatic_ring (0,1) $ do
  insertMolecule (atom C sp2) (0,1)
  newOrbital (6,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (6,2)
  insertMoleculeTwist (atom O sp3) (0,0) (-pi/6)
  hydrogenate

cyclohexane_chair = beginWith (atom C sp3) (0,0) $ do
  insertMoleculeTwist (atom C sp3) (0,0) (pi/3)
  newOrbital (1,1)
  insertMoleculeTwist (atom C sp3) (0,0) (pi/3)
  newOrbital (2,1)
  insertMoleculeTwist (atom C sp3) (0,0) (pi/3)
  newOrbital (3,2)
  insertMoleculeTwist (atom C sp3) (0,0) (pi/3)
  newOrbital (4,3)
  insertMoleculeTwist (atom C sp3) (0,0) (pi/3)
  matchEnds
  hydrogenate
  findConvexHull green

cyclohexane_boat = beginWith (atom C sp3) (0,0) $ do
  insertMoleculeTwist (atom C sp3) (0,0) (pi/3)
  newOrbital (1,1)
  insertMoleculeTwist (atom C sp3) (0,0) 0
  newOrbital (2,2)
  insertMoleculeTwist (atom C sp3) (0,0) (pi/3)
  newOrbital (3,1)
  insertMoleculeTwist (atom C sp3) (0,0) (pi/3)
  newOrbital (4,1)
  insertMoleculeTwist (atom C sp3) (0,0) 0
  matchEnds
  hydrogenate
  findConvexHull green

alpha_glucose_ring = beginWith (atom O sp3) (0,0) $ do
  insertMoleculeStretchTwist (atom C sp3) (0,0) 154 pi
  newOrbital (0,1)
  insertMoleculeStretchTwist (atom C sp3) (0,0) 154 pi
  newOrbital (1,3)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (2,1)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (3,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  matchEnds
  newOrbital (2,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (1,1)
  insertMoleculeTwist (atom O sp3) (0,0) (-2*pi/3)
  newOrbital (3,3)
  insertMoleculeTwist (atom O sp3) (0,0) 0
  newOrbital (4,2)
  insertMoleculeTwist (atom O sp3) (0,0) (-2*pi/3)
  newOrbital (5,1)
  insertMoleculeTwist (atom O sp3) (0,0) 0
  newOrbital (6,2)
  insertMoleculeTwist (atom O sp3) (0,0) pi
  hydrogenate
  findConvexHull green

beta_glucose_ring = beginWith (atom O sp3) (0,0) $ do
  insertMoleculeStretchTwist (atom C sp3) (0,0) 154 pi
  newOrbital (0,1)
  insertMoleculeStretchTwist (atom C sp3) (0,0) 154 pi
  newOrbital (1,3)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (2,1)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (3,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  matchEnds
  newOrbital (2,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (1,2)
  insertMoleculeTwist (atom O sp3) (0,0) (-2*pi/3)
  newOrbital (3,3)
  insertMoleculeTwist (atom O sp3) (0,0) 0
  newOrbital (5,1)
  insertMoleculeTwist (atom O sp3) (0,0) 0
  newOrbital (4,2)
  insertMoleculeTwist (atom O sp3) (0,0) (-2*pi/3)
  newOrbital (6,2)
  insertMoleculeTwist (atom O sp3) (0,0) pi
  hydrogenate
  findConvexHull green

edta = beginWith (atom C sp3) (0,0) $ do
  insertMolecule (atom C sp3) (0,0)
  newOrbital (0,1)
  insertMolecule (atom H s) (0,0)
  newOrbital (0,2)
  insertMolecule (atom H s) (0,0)
  newOrbital (1,1)
  insertMolecule (atom H s) (0,0)
  newOrbital (1,2)
  insertMolecule (atom H s) (0,0)
  newOrbital (0,3)

  insertMolecule (atom N sp3) (0,0)
  newOrbital (1,3)
  insertMolecule (atom N sp3) (0,0)

  newOrbital (6,1)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (8,2)
  insertMoleculeTwist (atom C sp2) (0,1) (-pi/2)
  newOrbital (9,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (9,2)
  insertMolecule (ion O (-1) sp3) (0,0)
  newOrbital (8,1)
  insertMolecule (atom H s) (0,0)
  newOrbital (8,3)
  insertMolecule (atom H s) (0,0)

  newOrbital (6,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (14,1)
  insertMoleculeTwist (atom C sp2) (0,1) (-pi/2)
  newOrbital (15,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (15,2)
  insertMolecule (ion O (-1) sp3) (0,0)
  newOrbital (14,2)
  insertMolecule (atom H s) (0,0)
  newOrbital (14,3)
  insertMolecule (atom H s) (0,0)

  newOrbital (7,1)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (20,2)
  insertMoleculeTwist (atom C sp2) (0,1) (-pi/2)
  newOrbital (21,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (21,2)
  insertMolecule (ion O (-1) sp3) (0,0)
  newOrbital (20,1)
  insertMolecule (atom H s) (0,0)
  newOrbital (20,3)
  insertMolecule (atom H s) (0,0)

  newOrbital (7,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (26,1)
  insertMoleculeTwist (atom C sp2) (0,1) (-pi/2)
  newOrbital (27,0)
  insertMolecule (atom O sp2) (0,0)
  newOrbital (27,2)
  insertMolecule (ion O (-1) sp3) (0,0)
  newOrbital (26,2)
  insertMolecule (atom H s) (0,0)
  newOrbital (26,3)
  insertMolecule (atom H s) (0,0)
  findConvexHull green

-- Useful Functional Groups

aromatic_ring = beginWith (atom C sp2) (0,0) $ do
  insertMoleculeTwist (atom C sp2) (0,0) pi
  newOrbital (1,2)
  insertMoleculeTwist (atom C sp2) (0,2) pi
  newOrbital (2,0)
  insertMoleculeTwist (atom C sp2) (0,0) pi
  newOrbital (3,2)
  insertMoleculeTwist (atom C sp2) (0,2) pi
  newOrbital (4,0)
  insertMoleculeTwist (atom C sp2) (0,0) pi
  matchEnds

phenyl = beginWith aromatic_ring (0,0) (hydrogenateExcept [(0,1)])

penta_C_sp3 = atomForceAngle C sp3 0 2 (3*pi/5)
penta_C_sp2 = atomForceAngle C sp2 0 2 (3*pi/5)

penta_hexa_C_sp2 = atomWarp C sp2 [ (0, Vector 0 (-71) 0)
                                  , (1, rotate3D (Vector 0 0 1) (3*pi/5) (Vector 0 (-71) 0))
                                  , (2, Vector (-sqrt 3 * 35.5) 35.5 0) ]

hexa_C_sp3 = atomForceAngle C sp3 0 3 (2*pi/3)

penta_N_sp3 = atomForceAngle N sp3 0 2 (3*pi/5)
penta_N_sp2 = atomForceAngle N sp2 0 2 (3*pi/5)

hexa_N_sp3 = atomForceAngle N sp3 0 2 (2*pi/3)

penta_O_sp3 = atomForceAngle O sp3 0 1 (3*pi/5)

hexa_O_sp3 = atomForceAngle O sp3 0 1 (2*pi/3)

-- Crystals

diamond = beginWith (atom C sp3) (0,0) $ do
  insertMoleculeTwist (atom C sp3) (0,0) pi
  let k = 308 / sqrt 3
  crystallize3D (Vector k k 0, 2)
                (Vector k 0 k, 3)
                (Vector 0 k k, 3)

                [ ( (1,0,0), [Bond (1,3) (0,3)] )
                , ( (0,1,0), [Bond (1,1) (0,2)] )
                , ( (0,0,1), [Bond (1,2) (0,1)] ) ]

_NaCl = beginWith (ion Na 1 emptyShell) (0,0) $ do
  let k = 564.02
  addMolecule (translate3D (Vector k 0 0) (ion Cl (-1) emptyShell))
  crystallize3D (Vector (2*k) 0 0, 5)
                (Vector    k  k 0, 5)
                (Vector    k  0 k, 5) []

cubic_NaCl = beginWith (ion Na 1 emptyShell) (0,0) $ do
  let k = 564.02
      l = 2 * k
  addMolecule (translate3D (Vector k 0 0) (ion Cl (-1) emptyShell))
  addMolecule (translate3D (Vector 0 k 0) (ion Cl (-1) emptyShell))
  addMolecule (translate3D (Vector 0 0 k) (ion Cl (-1) emptyShell))
  addMolecule (translate3D (Vector k k 0) (ion Na   1  emptyShell))
  addMolecule (translate3D (Vector 0 k k) (ion Na   1  emptyShell))
  addMolecule (translate3D (Vector k 0 k) (ion Na   1  emptyShell))
  addMolecule (translate3D (Vector k k k) (ion Cl (-1) emptyShell))
  crystallize3D (Vector l 0 0, 3)
                (Vector 0 l 0, 3)
                (Vector 0 0 l, 3) []

ice_Ih = beginWith (atom O sp3) (0,0) $ do
  let a = 451.81
      c = 735.60
      j = 2 * k
      k = a / 3
      l = c / sqrt 3
  insertMoleculeStretchTwist (atom O sp3) (0,2) (a / sqrt 3) (-pi/3)
  deleteBond (Bond (0,0) (1,2))
  newOrbital (0,3)
  insertMoleculeStretchTwist (atom O sp3) (0,1) (c / 2 - k / sqrt 3) (2*pi/3)
  deleteBond (Bond (0,3) (2,1))
  newOrbital (1,0)
  insertMoleculeStretchTwist (atom O sp3) (0,2) (c / 2 - k / sqrt 3) (2*pi/3)
  deleteBond (Bond (1,0) (3,2))
  hydrogenate
  addIMF 0 8
  addIMF 1 6
  addIMF 3 11
  crystallize3D (Vector (-l) (-l)   l, 2)
                (Vector   j    0    j, 3)
                (Vector   0    j    j, 3)

                [ ( (1,0,0), [IMF 2 10] )
                , ( (0,1,0) , [IMF 5 0] )
                , ( (0,0,1) , [IMF 1 4] )
                , ( (-1,1,0), [IMF 3 7] )
                , ( (-1,0,1), [IMF 9 2] ) ]

ice_Ic = beginWith (atom O sp3) (0,0) $ do
  let a = 635
      k = a / sqrt 2
  insertMoleculeStretchTwist (atom O sp3) (0,2) (a / sqrt 3) (-pi/3)
  deleteBond (Bond (0,0) (1,2))
  hydrogenate
  addIMF 5 1
  crystallize3D (Vector k k 0, 2)
                (Vector k 0 k, 3)
                (Vector 0 k k, 3)

                [ ( (1,0,0), [IMF 3 0] )
                , ( (0,1,0), [IMF 4 0] )
                , ( (0,0,1), [IMF 1 2] ) ]

graphene = beginWith (atom C sp2) (0,0) $ do
  insertMolecule (atom C sp2) (0,0)
  crystallize2D (Vector (-71 * sqrt 3) 213 0, 5)
                (Vector ( 71 * sqrt 3) 213 0, 5)
                [ ( (1,0), [Bond (0,2) (1,2)] )
                , ( (0,1), [Bond (0,1) (1,1)] ) ]

graphite = beginWith (atom C sp2) (0,0) $ do
  insertMolecule (atom C sp2) (0,0)
  addMolecule (translate3D (Vector (-71 * sqrt 3) 71 335) (atom C sp2))
  newOrbital (2,0)
  insertMolecule (atom C sp2) (0,0)
  crystallize3D (Vector (-71 * sqrt 3) 213 0, 5)
                (Vector ( 71 * sqrt 3) 213 0, 5)
                (Vector 0 0 670, 2)
                [ ( (1,0,0), [Bond (0,2) (1,2), Bond (2,2) (3,2), IMF 2 1] )
                , ( (0,1,0), [Bond (0,1) (1,1), Bond (2,1) (3,1)] )
                , ( (1,0,1), [IMF 2 1] ) ]

-- Polymers

bonding_alpha_glucose_ring = beginWith (atom O sp3) (0,0) $ do
  insertMoleculeStretchTwist (atom C sp3) (0,0) 154 pi
  newOrbital (0,1)
  insertMoleculeStretchTwist (atom C sp3) (0,0) 154 pi
  newOrbital (1,3)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (2,1)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (3,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  matchEnds
  newOrbital (2,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (1,1)
  insertMoleculeTwist (atomForceAngle O sp3 0 1 (2*pi/3 + pi/12)) (0,0) (-2*pi/3)
  newOrbital (3,3)
  insertMoleculeTwist (atom O sp3) (0,0) 0
  newOrbital (5,1)
  insertMoleculeTwist (atom O sp3) (0,0) 0
  newOrbital (6,2)
  insertMoleculeTwist (atom O sp3) (0,0) pi
  hydrogenateExcept [(4,2),(7,1)]
  findConvexHull green

amylose = beginWith bonding_alpha_glucose_ring (0,0) $ do
  polymerize (Bond (7,1) (4,2)) (-2*pi/3 - pi/9) [] 24

bonding_beta_glucose_ring = beginWith (atom O sp3) (0,0) $ do
  insertMoleculeStretchTwist (atom C sp3) (0,0) 154 pi
  newOrbital (0,1)
  insertMoleculeStretchTwist (atom C sp3) (0,0) 154 pi
  newOrbital (1,3)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (2,1)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (3,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  matchEnds
  newOrbital (2,2)
  insertMoleculeTwist (atom C sp3) (0,0) pi
  newOrbital (1,2)
  insertMoleculeTwist (atom O sp3) (0,0) (-2*pi/3)
  newOrbital (3,3)
  insertMoleculeTwist (atom O sp3) (0,0) 0
  newOrbital (5,1)
  insertMoleculeTwist (atom O sp3) (0,0) 0
  newOrbital (6,2)
  insertMoleculeTwist (atom O sp3) (0,0) pi
  hydrogenateExcept [(4,2),(7,1)]
  findConvexHull green

cellulose = beginWith bonding_beta_glucose_ring (0,0) $ do
  polymerize (Bond (7,1) (4,2)) (-2*pi/3) [] 12