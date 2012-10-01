
\section{Covariance Model Comparison}

\begin{code}
{-# LANGUAGE NamedFieldPuns, RecordWildCards, PatternGuards #-}
module Main where

-- Batteries
import qualified Data.MemoCombinators as Memo
import System.Environment (getArgs)
import Data.Array ((!),bounds)
import Debug.Trace
import Text.Regex.Posix
import Text.Printf
import Data.List
import Data.Ord
import Data.Function
import Data.Maybe

-- own Cabal
import Biobase.Infernal.CM
import Biobase.Infernal.CMImport
\end{code}



To be able to get the most information, do not use Double's but a datatype
holding more information.

\begin{code}
data Score = Score
  { score  :: Double -- total score
  , statel :: [Int] -- states traversed
  , stater :: [Int] -- states traversed
  , rnaseq :: String -- rna sequence "ACGU"
  } deriving (Show)

add :: Score -> Score -> Score
add sl sr = Score
  { score  = score sl + score sr
  , statel = statel sl ++ statel sr
  , stater = stater sl ++ stater sr
  , rnaseq = rnaseq sl ++ rnaseq sr
  }
\end{code}



Score calculation for two covariance models.

TODO: at each alternative, allow to use just the recursive case instead of
recursive+score to get a local model,string and score. That is:
max (rec+trans+emit) (rec)

\begin{code}
recurse rfl rfr local = rec 0 0 where
  --rec = Memo.memo2 Memo.integral Memo.integral rec'
  rec = Memo.memo2 (Memo.arrayRange (0,lastl)) (Memo.arrayRange (0,lastr)) rec'
  -- the recursive calculation, disallow loops into selfstates
  rec' :: Int -> Int -> Score
  rec' l r
    | tl == "E" && tr == "E" = Score 0 [sidl] [sidr] ""
    -- TODO arent we missing something here? like path scores?
    | tl == "E" = maxB [rec l (fst cr) | cr <- csr, r /= fst cr]
    | tr == "E" = maxB [rec (fst cl) r | cl <- csl, l /= fst cl]
    | tl == "B" && tr == "B" = maxB
      [ rec (left csl)  (left csr)  `add` rec (right csl) (right csr)
      , rec (left csl)  (right csr) `add` rec (right csl) lastr       `add` rec lastl (left csr)
      , rec (right csl) (left csr)  `add` rec (left csl)  lastr       `add` rec lastl (right csr)
      ]
    | tl == "B" = maxB
      [ rec (left csl)  r `add` rec (right csl) lastr
      , rec (right csl) r `add` rec (left csl)  lastr
      ]
    | tr == "B" = maxB
      [ rec l (left csr)  `add` rec lastl (right csr)
      , rec l (right csr) `add` rec lastl (left csr)
      ]
    | otherwise = maxB $ [rec cl cr | local, (cl,_) <- csl, (cr,_) <- csr, l /= cl || r /= cr] ++
                         [build cl cr | cl <- csl, cr <- csr, l /= fst cl || r /= fst cr]
    where
      tl = stype $ ssl ! l
      tr = stype $ ssr ! r
      csl = schildren $ ssl ! l
      csr = schildren $ ssr ! r
      esl = semissions $ ssl ! l
      esr = semissions $ ssr ! r
      sidl = sid $ ssl ! l
      sidr = sid $ ssr ! r
      infty = 1/0
      em4states = ["ML","IL","MR","IR"]
      sd = ["S","D"]
      -- extract transition score
      trans (_, Nothing) = 0
      trans (_, Just x) = x
      -- special lookup
      lookup' k vs = fromJust $ lookup k vs
      -- extract emission score
      emit
        | tl == tr && esl /= [] && esr /= [] =
          maximumBy (comparing snd) $ zipWith (\a b -> (fst a, snd a + snd b)) esl esr
        | esl /= [] && esr == [] = maximumBy (comparing snd) esl
        | esl == [] && esr /= [] = maximumBy (comparing snd) esr

        | tl /= tr && length esl == 16 && (tr == "IL" || tr == "ML") =
          maximumBy (comparing snd) $ map (\([b1,b2],sl) -> ([b1,b2],sl + lookup' [b1] esr)) esl
        | tl /= tr && length esr == 16 && (tl == "IL" || tl == "ML") =
          maximumBy (comparing snd) $ map (\([b1,b2],sr) -> ([b1,b2],sr + lookup' [b1] esl)) esr
        | tl /= tr && length esl == 16 && (tr == "IR" || tr == "MR") =
          maximumBy (comparing snd) $ map (\([b1,b2],sl) -> ([b1,b2],sl + lookup' [b2] esr)) esl
        | tl /= tr && length esr == 16 && (tl == "IR" || tl == "MR") =
          maximumBy (comparing snd) $ map (\([b1,b2],sr) -> ([b1,b2],sr + lookup' [b2] esl)) esr

        | tl `elem` em4states && tr `elem` em4states && drop 1 tl == drop 1 tr =
          maximumBy (comparing snd) $ zipWith (\a b -> (fst a, snd a + snd b)) esl esr
        | tl `elem` sd  && tr `elem` sd = ("",0)
        | drop 1 tl /= drop 1 tr =
            let
              (lb,ls) = maximumBy (comparing snd) esl
              (rb,rs) = maximumBy (comparing snd) esr
            in
              if drop 1 tl == "L"
              then (lb ++ rb, ls+rs)
              else (rb ++ lb, ls+rs)
        | otherwise = trace ("emit problem" ++ show (tl,tr)) ("",0)
      -- modify the rnasequence based on the given bases to add
      modifySeq s bs
        | [b1,b2] <- bs = [b1] ++ s ++ [b2]
        | tl == tr && tl `elem` sd = s
        | drop 1 tl == "L" = bs ++ s
        | drop 1 tl == "R" = s ++ bs
        | tl `elem` sd && drop 1 tr == "L" = bs ++ s
        | tl `elem` sd && drop 1 tr == "R" = s ++ bs
        | otherwise = error ("error in modifySeq: " ++ show (tl,tr,bs))
      -- build a complex step
      build cl cr = Score
                      { score  = score + trans cl + trans cr + es
                      , statel = [sidl] ++ statel
                      , stater = [sidr] ++ stater
                      , rnaseq = modifySeq rnaseq ebs
                      }
        where
          Score {..} = rec (fst cl) (fst cr)
          (ebs,es) = emit
  ssl = states rfl
  ssr = states rfr
  lastl = snd $ bounds ssl
  lastr = snd $ bounds ssr
  maxB = maximumBy (comparing score)
  left xs = fst $ xs !! 0
  right xs = fst $ xs !! 1
\end{code}



\begin{code}
main = do
  (a:b:[]) <- getArgs
  Right rf1 <- fromFile a
  Right rf2 <- fromFile b
  let only1 = (score $ recurse rf1 rf1) / 2
  let only2 = (score $ recurse rf2 rf2) / 2
  let reverse = recurse rf2 rf1
  let res = recurse rf1 rf2
  let relative = (score res - min only1 only2) / max only1 only2
  printf  "Relative: %6.2f Absolute: %7.2f %7.2f %7.2f Sequence: %s   Files: %s %s\n"
          relative
          (score res - min only1 only2)
          only1
          only2
          (rnaseq res)
          a
          b
\end{code}
