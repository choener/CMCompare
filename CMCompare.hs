{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | This program compares two Infernal covariance models with each other.
-- Based on the Infernal CM scoring mechanism, a Link sequence and Link score
-- are calculated. The Link sequence is defined as the sequence scoring highest
-- in both models simultanuously.
--
-- The complete algorithm is described in:
--
-- "Christian Höner zu Siederdissen, and Ivo L. Hofacker. 2010. Discriminatory
-- power of RNA family models. Bioinformatics 26, no. 18: 453–59."
--
-- <http://bioinformatics.oxfordjournals.org/content/26/18/i453.long>
--
--
--
-- NOTE always use coverage analysis to find out, if we really used all code
-- paths (in long models, if a path is not taken, there is a bug)

-- NOTE when comparing hits with cmsearch, use the following commandline:
--
-- cmsearch --no-null3 --cyk --fil-no-hmm --fil-no-qdb
--
-- --no-null3 : important, the test sequence is so short that null3 can easily
-- generate scores that are way off! remember, we are interested in a sequence
-- that is typically embedded in something large
--
-- --fil-no-hmm, --fil-no-qdb: do not use heuristics for speedup, they
-- sometimes hide results (in at least one case)
--
-- (--toponly): if the comparison was done onesided
--
-- (-g): if you want to compare globally



module Main where

import Control.Arrow (first,second,(***))
import Control.Monad
import Data.Array.IArray
import System.Console.CmdArgs
import Text.Printf

import Biobase.SElab.CM
import Biobase.SElab.CM.Import
import Biobase.SElab.Types

import BioInf.CMCompare



-- * Handling user I/O.
--
-- TODO add an option to filter by minimal score (default: -1000000) if the
-- filter is on, we get a result only, if the score is above the threshold

data Options = Options
  { global :: Bool
  , invert :: Bool
  , pbegin :: Double
  , pend :: Double
  , endmarker :: Bool
  , nobeginilir :: Bool
  , fastIns :: Bool
  , models :: [String]
  } deriving (Show,Data,Typeable)

options = Options
  { global      = False &= help "global model comparison"
  , invert      = False &= help "TODO invert second model (read right to left)"
  , pbegin      = 0.05  &= help "aggregate local begin probability"
  , pend        = 0.05  &= help "aggregate local end probability"
  , endmarker   = False &= help "add an endmarker into the rnastring to denote local ends"
  , nobeginilir = False &= help "trailing left or right nucleotides change the score"
  , fastIns     = False &= help "fast insertion heuristic"
  , models      = def   &= args -- &= help "path to exactly two covariance models"
  } &= summary "CMCompare: Discriminatory Power of RNA Family Models" &= help "(c) 2010, Christian Hoener zu Siederdissen and Ivo Hofacker\nchoener@tbi.univie.ac.at\nLicensed under the GPLv3\n" &= verbosity
-- TODO put fixTransition in BiobaseInfernal, rename and whatnot

-- | what is this? this sets the transition score in the root to 0 for the
-- transition into IL and IR. For CM comparison, this makes sense. Consider two
-- cms .(.) and (.); given "acag" both should score about the same, but the
-- second can only, if IL eats a nucleotide. for single CM search, this is
-- taken care of by the DP algorithm which doesn't work here.

{-
fixTransition :: CM' -> CM'
fixTransition x = x{states = ss // [(0,rtnew)]} where
  ss = states x
  rt = ss ! 0
  tr = schildren rt
  trnew = [Transition 1 0, Transition 2 0] ++ drop 2 tr
  rtnew = rt{schildren = trnew}
-}

applyIf c f = if c then f else id

answers fastIns optf m1 m2 = map (finalize *** finalize) . opt . concat . elems $ recurse fastIns optf m1 m2 where
  (end,lbegin,start,delete,matchP,matchL,insertL,matchR,insertR,branch,opt,finalize) = optf

results fastIns optf m1 m2 = opt . concat . elems $ recurse fastIns optf m1 m2 where
  (end,lbegin,start,delete,matchP,matchL,insertL,matchR,insertR,branch,opt,finalize) = optf

main = do
  Options{..} <- cmdArgs options
  normal <- isNormal
  loud <- isLoud
  let quiet = not normal
  unless (length models == 2) $ do
    fail "give exactly two CMs"
  let [a,b] = models
  [theA] <- fromFile a
  [theB] <- fromFile b
  let m1 = if nobeginilir then theA else makeLocal pbegin pend theA -- TODO !!! applyIf (not nobeginilir) fixTransition $ applyIf (not global) (cmMakeLocal pbegin pend) m1'
  let m2 = if nobeginilir then theA else makeLocal pbegin pend theB -- TODO !!! applyIf (not nobeginilir) fixTransition $ applyIf (not global) (cmMakeLocal pbegin pend) m2'
  let pr = (\(x,y) -> putStrLn x >> putStrLn "" >> putStrLn y >> putStrLn "++++++++++++")
  when (quiet && not normal) $ do
    let (a1,a2) = head $ results fastIns cykMaxiMin m1 m2
    printf "%s   %s %10.3f %10.3f\n" a b (unBitScore a1) (unBitScore a2)
  when (normal && not loud) $ do
    let ((((cyk1,vn1),rna1),db1),(((cyk2,vn2),_),db2)) = head $ results fastIns (cykMaxiMin <*> visitedNodes <*> rnaString endmarker <*> dotBracket endmarker) m1 m2
    let rnaBoth = map f rna1 where
          f x
            | x=='N' = '_'
            | otherwise   = x
    printf "%s   %s %10.3f %10.3f %s %s %s %s %s\n" a b (unBitScore cyk1) (unBitScore cyk2) rnaBoth db1 db2 (show $ map unNodeID vn1) (show $ map unNodeID vn2)
  when loud . mapM_ pr $ answers fastIns (cykMaxiMin <*> visitedNodes <*> rnaString endmarker <*> dotBracket endmarker <*> extendedOutput) m1 m2

