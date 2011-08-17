module Algorithms where

import Geometry
import Chemistry

import Data.List

-- Warping Atoms

bondRepulsion :: [Vector] -> [Vector] -> Int -> [Vector]
bondRepulsion _ [] _ = []
bondRepulsion _ free 0 = free
bondRepulsion fixed free n = bondRepulsion fixed (map (repel (fixed ++ free)) free) (n - 1)
  where repel vs v = (v + sum forces) `scaleTo` (norm v)
          where forces = map (repulsion v) (delete v vs)
                repulsion v u = v - u `scaleTo` (100 / (angleBetween v u)^2)

warpAtom :: Atom -> [(Int, Vector)] -> Atom
warpAtom (Atom p e q h bs) ivs =
  let (ibs1, ibs2, fixed, free) = separate (zip [0..] bs) ivs
      ibs1' = zipWith (\(i, BondSite os bt u rx) v -> (i, BondSite os bt v (alignTo v u rx))) ibs1 fixed
      ibs2' = zipWith (\(i, BondSite os bt u rx) v -> (i, BondSite os bt v (alignTo v u rx))) ibs2 (bondRepulsion fixed free 1000)
  in Atom p e q h (regroup ibs1' ibs2')
  where separate (ib@(i, b):ibs) ivs =
          let (ibs1, ibs2, fixed, free) = separate ibs ivs
              v = direction b
          in case lookup i ivs of
               Just v -> (ib:ibs1,    ibs2, v:fixed,   free)
               _      -> (   ibs1, ib:ibs2,   fixed, v:free)
        separate _ _ = ([],[],[],[])

        regroup xs@((i1,b1):ibs1) ys@((i2,b2):ibs2) =
          if i1 < i2
           then b1 : regroup ibs1 ys
           else b2 : regroup xs ibs2
        regroup ibs1 [] = map snd ibs1
        regroup [] ibs2 = map snd ibs2
