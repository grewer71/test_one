{-# LANGUAGE RecordWildCards #-}

import Numeric.LinearAlgebra (Matrix, R, (<>), inv, ident, sym, toRows, fromLists)
import Numeric.LinearAlgebra.Data (loadMatrix)
import Numeric.LinearAlgebra.Algorithms (eigSH)
import System.IO (readFile)
import Data.List (transpose)

-- Решение уравнения Риккати: AᵀP + PA - PBR⁻¹BᵀP + Q = 0
solveRiccati :: Matrix R -> Matrix R -> Matrix R -> Matrix R -> Matrix R
solveRiccati a b q r = p
  where
    -- Размерности
    n = rows a
    m = cols b
    
    -- Вспомогательные матрицы
    rInv = inv r
    bTrans = tr b
    aTrans = tr a
    
    -- Формируем гамильтонову матрицу
    ham = fromBlocks [[a, -b <> rInv <> bTrans], 
                      [-q, -aTrans]]
    
    -- Вычисляем собственные векторы
    (_, vecs) = eigSH (sym ham)
    
    -- Берем первые n собственных векторов
    u = takeRows n vecs
    v = dropRows n vecs
    
    -- Решение P = V * U⁻¹
    p = v <> inv u

-- Вспомогательные функции
fromBlocks :: [[Matrix R]] -> Matrix R
fromBlocks = fromLists . concatMap (map toList . toRows) . concat . transpose

takeRows :: Int -> Matrix R -> Matrix R
takeRows k = fromLists . take k . toLists

dropRows :: Int -> Matrix R -> Matrix R
dropRows k = fromLists . drop k . toLists

main :: IO ()
main = do
    -- Чтение матриц из файлов
    a <- loadMatrix "a.txt"
    b <- loadMatrix "b.txt"
    q <- loadMatrix "q.txt"  -- или "e.txt"
    
    -- Матрица R (единичная размерности cols b)
    let r = ident (cols b)
    
    -- Решение уравнения Риккати
    let p = solveRiccati a b q r
    
    -- Вывод результата
    putStrLn "Решение P:"
    print p
    
    -- Проверка невязки
    let residual = tr a <> p + p <> a - p <> b <> inv r <> tr b <> p + q
    putStrLn "\nНевязка (должна быть близка к нулю):"
    print residual